package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**
 * Assembles fused rings named using fusion nomenclature for outputParse.
 * This is a separate implementation that does not detach/remove groups during processing.
 * @author dl387
 *
 */
class FusedRingBuilderForOutput {
	private final BuildState state;
	private final List<Element> groupsInFusedRing;
	private final Element lastGroup;
	private final Fragment parentRing;
	private final Map<Integer,Fragment> fragmentInScopeForEachFusionLevel = new HashMap<>();
	private final Map<Atom, Atom> atomsToRemoveToReplacementAtom = new HashMap<>();
	// MolLangData: Map to track original labels for each atom before fusion
	// Maps: Atom -> (ring index in groupsInFusedRing, original label from that ring)
	private final Map<Atom, Map<Integer, String>> atomToOriginalLabels = new HashMap<>();

	private FusedRingBuilderForOutput(BuildState state, List<Element> groupsInFusedRing) {
		this.state = state;
		this.groupsInFusedRing = groupsInFusedRing;
		lastGroup = groupsInFusedRing.get(groupsInFusedRing.size()-1);
		parentRing = lastGroup.getFrag();
		fragmentInScopeForEachFusionLevel.put(0, parentRing);
	}

	/**
	 * Master method for processing fused rings. Fuses groups together
	 * @param state: contains the current id and fragment manager
	 * @param subOrRoot Element (substituent or root)
	 * @throws StructureBuildingException
	 */
	static void processFusedRings(BuildState state, Element subOrRoot) throws  StructureBuildingException {
		List<Element> groups = subOrRoot.getChildElements(GROUP_EL);
		if (groups.size() < 2){
			return;//nothing to fuse
		}
		List<Element> groupsInFusedRing =new ArrayList<>();
		for (int i = groups.size()-1; i >=0; i--) {//group groups into fused rings
			Element group =groups.get(i);
			groupsInFusedRing.add(0, group);
			if (i!=0){
				Element startingEl = group;
				if ((group.getValue().equals("benz") || group.getValue().equals("benzo")) && FUSIONRING_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
					Element beforeBenzo = OpsinTools.getPreviousSibling(group);
					if (beforeBenzo !=null && beforeBenzo.getName().equals(LOCANT_EL)){
						startingEl = beforeBenzo;
					}
				}
				Element possibleGroup = OpsinTools.getPreviousSiblingIgnoringCertainElements(startingEl, new String[]{MULTIPLIER_EL, FUSION_EL});
				if (!groups.get(i-1).equals(possibleGroup)){//end of fused ring system
					if (groupsInFusedRing.size()>=2){
						//This will be invoked in cases where there are multiple fused ring systems in the same subOrRoot such as some spiro systems
						new FusedRingBuilderForOutput(state, groupsInFusedRing).buildFusedRing();
					}
					groupsInFusedRing.clear();
				}
			}
		}
		if (groupsInFusedRing.size()>=2){
			new FusedRingBuilderForOutput(state, groupsInFusedRing).buildFusedRing();
		}
	}

	/**
	 * Combines the groups given in the {@link FusedRingBuilderForOutput} constructor to destructively create the fused ring system
	 * This fused ring is then numbered
	 * @throws StructureBuildingException
	 */
	void buildFusedRing() throws StructureBuildingException{
		/*
		 * Apply any nonstandard ring numbering, sorts atomOrder by locant
		 * Aromatises appropriate cycloalkane rings, Rejects groups with acyclic atoms
		 */
        processRingNumberingAndIrregularities();
		// MolLangData: Initialize original labels mapping before fusion
		//initializeOriginalLabelsMapping();
		processBenzoFusions();//FR-2.2.8  e.g. in 2H-[1,3]benzodioxino[6',5',4':10,5,6]anthra[2,3-b]azepine  benzodioxino is one component
		List<Element> nameComponents = formNameComponentList();
		nameComponents.remove(lastGroup);

		List<Fragment> componentFragments = new ArrayList<>();//all the ring fragments (other than the parentRing). These will later be merged into the parentRing
		List<Fragment> parentFragments = new ArrayList<>();
		parentFragments.add(parentRing);
		
		int numberOfParents = 1;
		Element possibleMultiplier = OpsinTools.getPreviousSibling(lastGroup);
		if (nameComponents.size()>0 && possibleMultiplier !=null && possibleMultiplier.getName().equals(MULTIPLIER_EL)){
			numberOfParents = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
			possibleMultiplier.detach();
			for (int j = 1; j < numberOfParents; j++) {
				Fragment copyOfParentRing =state.fragManager.copyFragment(parentRing);
				parentFragments.add(copyOfParentRing);
				componentFragments.add(copyOfParentRing);
			}
		}

		/*The indice from nameComponents to use next. Work from right to left i.e. starts at nameComponents.size()-1*/
		int ncIndice = processMultiParentSystem(parentFragments, nameComponents, componentFragments);//handle multiparent systems
		/*
		 * The number of primes on the component to be connected. 
		 * This is initially 0 indicating fusion of unprimed locants with the letter locants of the parentRing
		 * Subsequently it will switch to 1 indicating fusion of a second order component (primed locants) with a 
		 * first order component (unprimed locants)
		 * Next would be double primed fusing to single primed locants etc.
		 * 
		 */
		int fusionLevel = (nameComponents.size()-1 -ncIndice)/2;
		for (; ncIndice>=0; ncIndice--) {
			Element fusion = null;
			if (nameComponents.get(ncIndice).getName().equals(FUSION_EL)){
				fusion = nameComponents.get(ncIndice--);
			}
			if (ncIndice <0 || !nameComponents.get(ncIndice).getName().equals(GROUP_EL)){
				throw new StructureBuildingException("Group not found where group expected. This is probably a bug");
			}
			Fragment nextComponent = nameComponents.get(ncIndice).getFrag();
			int multiplier = 1;
			Element possibleMultiplierEl = OpsinTools.getPreviousSibling(nameComponents.get(ncIndice));//e.g. the di of difuro
			if (possibleMultiplierEl != null && possibleMultiplierEl.getName().equals(MULTIPLIER_EL)){
				multiplier = Integer.parseInt(possibleMultiplierEl.getAttributeValue(VALUE_ATR));
			}
			String[] fusionDescriptors =null;
			if (fusion !=null){
				String fusionDescriptorString = fusion.getValue().toLowerCase(Locale.ROOT).substring(1, fusion.getValue().length()-1);
				if (multiplier ==1){
					fusionDescriptors = new String[]{fusionDescriptorString};
				}
				else{
					if (fusionDescriptorString.split(";").length >1){
						fusionDescriptors = fusionDescriptorString.split(";");
					}
					else if (fusionDescriptorString.split(":").length >1){
						fusionDescriptors = fusionDescriptorString.split(":");
					}
					else if (fusionDescriptorString.split(",").length >1){
						fusionDescriptors = fusionDescriptorString.split(",");
					}
					else{//multiplier does not appear to mean multiplied component. Could be indicating multiplication of the whole fused ring system
						if (ncIndice!=0){
							throw new StructureBuildingException("Unexpected multiplier: " + possibleMultiplierEl.getValue() +" or incorrect fusion descriptor: " + fusionDescriptorString);
						}
						multiplier =1;
						fusionDescriptors = new String[]{fusionDescriptorString};
					}
				}
			}
			if (multiplier >1){
				possibleMultiplierEl.detach();
			}
			Fragment[] fusionComponents = new Fragment[multiplier];
			for (int j = 0; j < multiplier; j++) {
				if (j>0){
					fusionComponents[j] = state.fragManager.copyAndRelabelFragment(nextComponent,  j);
				}
				else{
					fusionComponents[j] = nextComponent;
				}
			}
			
			for (int j = 0; j < multiplier; j++) {
				Fragment component = fusionComponents[j];
				componentFragments.add(component);
				if (fusion !=null){
					if (fusionDescriptors[j].split(":").length==1){//A fusion bracket without a colon is used when applying to the parent component (except in a special case where locants are ommitted)
						//check for case of omitted locant from a higher order fusion bracket e.g. cyclopenta[4,5]pyrrolo[2,3-c]pyridine
						if (fusionDescriptors[j].split("-").length==1 && 
								fusionDescriptors[j].split(",").length >1 &&
								FragmentTools.allAtomsInRingAreIdentical(component)
								&& ((StringTools.countTerminalPrimes(fusionDescriptors[j].split(",")[0])) != fusionLevel) ){//Could be like cyclopenta[3,4]cyclobuta[1,2]benzene where the first fusion to occur has parent locants omitted not child locants
							int numberOfPrimes = StringTools.countTerminalPrimes(fusionDescriptors[j].split(",")[0]);
							//note that this is the number of primes on the parent ring. So would expect the child ring and hence the fusionLevel to be 1 higher
							if (numberOfPrimes + 1 != fusionLevel){
								if (numberOfPrimes + 2 == fusionLevel){//ring could be in previous fusion level e.g. the benzo in benzo[10,11]phenanthro[2',3',4',5',6':4,5,6,7]chryseno[1,2,3-bc]coronene
									fusionLevel--;
								}
								else{
									throw new StructureBuildingException("Incorrect number of primes in fusion bracket: " +fusionDescriptors[j]);
								}
							}
							relabelAccordingToFusionLevel(component, fusionLevel);
							List<String> numericalLocantsOfParent = Arrays.asList(fusionDescriptors[j].split(","));
							List<String> numericalLocantsOfChild = findPossibleNumericalLocants(component, determineAtomsToFuse(fragmentInScopeForEachFusionLevel.get(fusionLevel), numericalLocantsOfParent, null).size()-1);
							processHigherOrderFusionDescriptors(component, fragmentInScopeForEachFusionLevel.get(fusionLevel), numericalLocantsOfChild, numericalLocantsOfParent);
						}
						else{
							fusionLevel = 0;
							relabelAccordingToFusionLevel(component, fusionLevel);
							String fusionDescriptor = fusionDescriptors[j];
							String[] fusionArray = determineNumericalAndLetterComponents(fusionDescriptor);
							int numberOfPrimes =0;
							if (!fusionArray[1].equals("")){
								numberOfPrimes =StringTools.countTerminalPrimes(fusionArray[1]);
								if (fusionArray[0].equals("")){
									fusionDescriptor = fusionArray[1].replaceAll("'", "");
								}
								else{
									fusionDescriptor = fusionArray[0]+ "-" +fusionArray[1].replaceAll("'", "");
								}
								if (numberOfPrimes >= parentFragments.size()){
									throw new StructureBuildingException("Unexpected prime in fusion descriptor");
								}
							}
							performSimpleFusion(fusionDescriptor, component, parentFragments.get(numberOfPrimes));//e.g. pyrano[3,2-b]imidazo[4,5-e]pyridine where both are level 0 fusions
						}
					}
					else{
						//determine number of primes in fusor and hence determine fusion level
						int numberOfPrimes = -j + StringTools.countTerminalPrimes(fusionDescriptors[j].split(",")[0]);
						if (numberOfPrimes != fusionLevel){
							if (fusionLevel == numberOfPrimes +1){
								fusionLevel--;
							}
							else{
								throw new StructureBuildingException("Incorrect number of primes in fusion bracket: " +fusionDescriptors[j]);
							}
						}
						relabelAccordingToFusionLevel(component, fusionLevel);
						performHigherOrderFusion(fusionDescriptors[j], component, fragmentInScopeForEachFusionLevel.get(fusionLevel));
					}
				}
				else{
					relabelAccordingToFusionLevel(component, fusionLevel);
					performSimpleFusion(null, component, fragmentInScopeForEachFusionLevel.get(fusionLevel));
				}
			}
			fusionLevel++;
			if (multiplier ==1){//multiplied components may not be substituted onto
				fragmentInScopeForEachFusionLevel.put(fusionLevel, fusionComponents[0]);
			}
		}

		// MolLangData: Create snapshots of labels and atoms before fusion
		//Map<Atom, String> parentRingAtomToLabel = createAtomToLabelSnapshot(parentRing);
		// Have a list of mapping of all the componentsFragment to their atomToLabel snapshot
		//List<Map<Atom, String>> componentFragmentsAtomToLabel = new ArrayList<>();

		for (Fragment ring : componentFragments) {
			//Map<Atom, String> componentFragmentAtomToLabel = createAtomToLabelSnapshot(ring);
			state.fragManager.incorporateFragment(ring, parentRing);
			//updateSnapshotWithAtomReplacements(parentRingAtomToLabel, atomsToRemoveToReplacementAtom);
			//updateSnapshotWithAtomReplacements(componentFragmentAtomToLabel, atomsToRemoveToReplacementAtom);
			//componentFragmentsAtomToLabel.add(componentFragmentAtomToLabel);
		}
		
		// MolLangData: Update atomToOriginalLabels to map removed atoms to their replacement atoms
		//updateAtomToOriginalLabelsWithReplacements(atomsToRemoveToReplacementAtom);

		
		removeMergedAtoms();

		
		
		// MolLangData: Update snapshots to map removed atoms to their replacement atoms
		//updateSnapshotWithAtomReplacements(parentRingAtomToLabel, atomsToRemoveToReplacementAtom);
		//updateSnapshotWithAtomReplacements(fusedRingAtomToLabel, atomsToRemoveToReplacementAtom);

		FusedRingNumberer.numberFusedRing(parentRing);//numbers the fused ring;
		
		// MolLangData: create a fusedRingNumbering element to store the numbering information
		//createFusedRingNumberingElement(parentRing, parentRingAtomToLabel, componentFragmentsAtomToLabel);

		StringBuilder fusedRingName = new StringBuilder();
		for (Element element : nameComponents) {
			fusedRingName.append(element.getValue());
		}
		fusedRingName.append(lastGroup.getValue());

		Element fusedRingEl =lastGroup;//reuse this element to save having to remap suffixes...

		/*
		// MolLangData: we copy the original tokenEl to a new fusedChildRing element, and add it at the first position of the parent ring tokenEl
		List<Attribute> allAttributes = fusedRingEl.getAttributes();
		// Create a new fusedChildRing element
		Element fusedChildRingEl = new TokenEl(FUSEDCHILDRING_EL);
		// Set the value as the copy of the original tokenEl value
		fusedChildRingEl.setValue(fusedRingEl.getValue());
		// Copy all attributes from the original tokenEl to the new fusedChildRing element
		for (Attribute attribute : allAttributes) {
			fusedChildRingEl.addAttribute(new Attribute(attribute));
		}
		// insert the new fusedChildRing element at the first position of the parent ring tokenEl
		fusedRingEl.insertChild(fusedChildRingEl, 0);
		*/

		fusedRingEl.getAttribute(VALUE_ATR).setValue(fusedRingName.toString()); //MolLangData Comment: this is original code, set the value of the fused ring as the original tokenEl value

		/* 
		List<Attribute> attributesToRemove = new ArrayList<>();
		for (Attribute attribute : allAttributes) {
			if (!attribute.getName().equals(VALUE_ATR) && !attribute.getName().equals(TYPE_ATR) && !attribute.getName().equals(SUBTYPE_ATR)) {
				attributesToRemove.add(attribute);
			}
		}
		for (Attribute attribute : attributesToRemove) {
			fusedRingEl.removeAttribute(attribute);
		}
		*/
		fusedRingEl.getAttribute(TYPE_ATR).setValue(RING_TYPE_VAL);
		fusedRingEl.setValue(fusedRingName.toString());

		for (Element element : nameComponents) {
			element.detach();
		}
	}

	private void removeMergedAtoms() {
		for (Atom a : atomsToRemoveToReplacementAtom.keySet()) {
			state.fragManager.removeAtomAndAssociatedBonds(a);
		}
		atomsToRemoveToReplacementAtom.clear();
	}

	/**
	 * Forms a list a list of all group and fusion elements between the first and last group in the fused ring
	 * @return
	 */
	private List<Element> formNameComponentList() {
		List<Element> nameComponents  = new ArrayList<>();
		Element currentEl = groupsInFusedRing.get(0);
		while(currentEl != lastGroup){
			if (currentEl.getName().equals(GROUP_EL) || currentEl.getName().equals(FUSION_EL)){
				nameComponents.add(currentEl);
			}
			currentEl = OpsinTools.getNextSibling(currentEl);
		}
		return nameComponents;
	}

	private void processRingNumberingAndIrregularities() throws StructureBuildingException {
		for (Element group : groupsInFusedRing) {
            Fragment ring = group.getFrag();
            if (ALKANESTEM_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
            	aromatiseCyclicAlkane(group);
            }
            processPartiallyUnsaturatedHWSystems(group, ring);
            if (group == lastGroup) {
                //perform a quick check that every atom in this group is infact cyclic. Fusion components are enumerated and hence all guaranteed to be purely cyclic
                List<Atom> atomList = ring.getAtomList();
                for (Atom atom : atomList) {
                    if (!atom.getAtomIsInACycle()) {
                        throw new StructureBuildingException("Inappropriate group used in fusion nomenclature. Only groups composed entirely of atoms in cycles may be used. i.e. not: " + group.getValue());
                    }
                }
                if (group.getAttribute(FUSEDRINGNUMBERING_ATR) != null) {
					// MolLangData: have a debug warning here, we have no idea what is this for, we should take care when this is happening
					state.addWarning(OpsinWarning.OpsinWarningType.MolLangData_DEBUG_WARNING, "MolLangData does not know when the code can run into this situation for this fused ring numbering; please report this issue to the developers");

                    String[] standardNumbering = group.getAttributeValue(FUSEDRINGNUMBERING_ATR).split("/", -1);
                    for (int j = 0; j < standardNumbering.length; j++) {
                        atomList.get(j).replaceLocants(standardNumbering[j]);
                    }
                } else {
                    ring.sortAtomListByLocant();//for those where the order the locants are in is sensible					}
                }
                for (Atom atom : atomList) {
                    atom.clearLocants();//the parentRing does not have locants, letters are used to indicate the edges
                }
            } else if (group.getAttribute(FUSEDRINGNUMBERING_ATR) == null) {
                ring.sortAtomListByLocant();//for those where the order the locants are in is sensible
            }
        }
	}

	/**
	 * Interprets the unlocanted unsaturator after a partially unsaturated HW Rings as indication of spare valency and detaches it
     * This is necessary as this unsaturator can only refer to the HW ring and for names like 2-Benzoxazolinone to avoid confusion as to what the 2 refers to.
	 * @param group
	 * @param ring
	 */
	private void processPartiallyUnsaturatedHWSystems(Element group, Fragment ring) {
		if (HANTZSCHWIDMAN_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR)) && group.getAttribute(ADDBOND_ATR)!=null){
			List<Element> unsaturators = OpsinTools.getNextAdjacentSiblingsOfType(group, UNSATURATOR_EL);
			if (unsaturators.size()>0){
				Element unsaturator = unsaturators.get(0);
				if (unsaturator.getAttribute(LOCANT_ATR)==null && unsaturator.getAttributeValue(VALUE_ATR).equals("2")){
					// MolLangData: have a debug warning here, we have no idea what is this for, we should take care when this is happening
					state.addWarning(OpsinWarning.OpsinWarningType.MolLangData_DEBUG_WARNING, "MolLangData does not know when the code can run into this situation for this unsaturator; please report this issue to the developers");

					unsaturator.detach();
					List<Bond> bondsToUnsaturate = StructureBuildingMethods.findBondsToUnSaturate(ring, 2, true);
					if (bondsToUnsaturate.isEmpty()) {
						throw new RuntimeException("Failed to find bond to unsaturate on partially saturated HW ring");
					}
					Bond b = bondsToUnsaturate.get(0);
					b.getFromAtom().setSpareValency(true);
					b.getToAtom().setSpareValency(true);
				}
			}
		}		
	}

	/**
	 * Given a cyclicAlkaneGroup determines whether or not it should be aromatised. Unlocanted ene will be detached if it is an aromatisation hint
	 * No unsaturators -->aromatise
	 * Just ane -->don't
	 * More than 1 ene or locants on ene -->don't
	 * yne --> don't
	 * @param cyclicAlkaneGroup
	 */
	private void aromatiseCyclicAlkane(Element cyclicAlkaneGroup) {
		Element next = OpsinTools.getNextSibling(cyclicAlkaneGroup);
		List<Element> unsaturators = new ArrayList<>();
		while (next!=null && next.getName().equals(UNSATURATOR_EL)){
			// MolLangData: we should have a debug warning here, we should take care when this is happening
			state.addWarning(OpsinWarning.OpsinWarningType.MolLangData_DEBUG_WARNING, "MolLangData does not know when the code can run into this situation for this unsaturator; please report this issue to the developers");
			unsaturators.add(next);
			next = OpsinTools.getNextSibling(next);
		}
		boolean conjugate =true;
		if (unsaturators.size()==1){
			int value = Integer.parseInt(unsaturators.get(0).getAttributeValue(VALUE_ATR));
			if (value !=2){
				conjugate =false;
			}
			else if (unsaturators.get(0).getAttribute(LOCANT_ATR)!=null){
				conjugate =false;
			}
		}
		else if (unsaturators.size()==2){
			int value1 = Integer.parseInt(unsaturators.get(0).getAttributeValue(VALUE_ATR));
			if (value1 !=1){
				conjugate =false;
			}
			else{
				int value2 = Integer.parseInt(unsaturators.get(1).getAttributeValue(VALUE_ATR));
				if (value2 !=2 || unsaturators.get(1).getAttribute(LOCANT_ATR)!=null){
					conjugate =false;
				}
			}
		}
		else if (unsaturators.size() >2){
			conjugate =false;
		}
		if (conjugate){
			for (Element unsaturator : unsaturators) {
				unsaturator.detach();
			}
			List<Atom> atomList = cyclicAlkaneGroup.getFrag().getAtomList();
		    for (Atom atom : atomList) {
		        atom.setSpareValency(true);
		    }
		}
	}

	private int processMultiParentSystem(List<Fragment> parentFragments, List<Element> nameComponents, List<Fragment> componentFragments) throws StructureBuildingException {
		int i = nameComponents.size()-1;
		int fusionLevel =0;
		if (i>=0 && parentFragments.size()>1){
			List<Fragment> previousFusionLevelFragments = parentFragments;
			for (; i>=0; i--) {
				if (previousFusionLevelFragments.size()==1){//completed multi parent system
					fragmentInScopeForEachFusionLevel.put(fusionLevel, previousFusionLevelFragments.get(0));
					break;
				}
				Element fusion = null;
				if (nameComponents.get(i).getName().equals(FUSION_EL)){
					fusion = nameComponents.get(i--);
				}
				else{
					throw new StructureBuildingException("Fusion bracket not found where fusion bracket expected");
				}
				if (i <0 || !nameComponents.get(i).getName().equals(GROUP_EL)){
					throw new StructureBuildingException("Group not found where group expected. This is probably a bug");
				}
				Fragment nextComponent = nameComponents.get(i).getFrag();
				relabelAccordingToFusionLevel(nextComponent, fusionLevel);
				int multiplier = 1;
				Element possibleMultiplierEl = OpsinTools.getPreviousSibling(nameComponents.get(i));
				if (possibleMultiplierEl != null && possibleMultiplierEl.getName().equals(MULTIPLIER_EL)){
					// MolLangData: we should have a debug warning here, we should take care when this is happening
					state.addWarning(OpsinWarning.OpsinWarningType.MolLangData_DEBUG_WARNING, "MolLangData does not know when the code can run into this situation for this multiplier; please report this issue to the developers");
					multiplier = Integer.parseInt(possibleMultiplierEl.getAttributeValue(VALUE_ATR));
					possibleMultiplierEl.detach();
				}
				List<Fragment> fusionComponents = new ArrayList<>();
				for (int j = 0; j < multiplier; j++) {
					if (j>0){
						Fragment clonedFrag = state.fragManager.copyFragment(nextComponent);
						relabelAccordingToFusionLevel(clonedFrag, j);//fusionLevels worth of primes already added
						fusionComponents.add(clonedFrag);
					}
					else{
						fusionComponents.add(nextComponent);
					}
				}
				fusionLevel+=multiplier;
				if (multiplier>1 && multiplier != previousFusionLevelFragments.size()){
					throw new StructureBuildingException("Mismatch between number of components and number of parents in fused ring system");
				}
				String fusionDescriptorString = fusion.getValue().toLowerCase(Locale.ROOT).substring(1, fusion.getValue().length()-1);
				String[] fusionDescriptors =null;
				if (fusionDescriptorString.split(";").length >1){
					fusionDescriptors = fusionDescriptorString.split(";");
				}
				else if (fusionDescriptorString.split(":").length >1){
					fusionDescriptors = fusionDescriptorString.split(":");
				}
				else if (fusionDescriptorString.split(",").length >1){
					fusionDescriptors = fusionDescriptorString.split(",");
				}
				else{
					throw new StructureBuildingException("Invalid fusion descriptor: " + fusionDescriptorString);
				}
				if (fusionDescriptors.length != previousFusionLevelFragments.size()){
					throw new StructureBuildingException("Invalid fusion descriptor: "+fusionDescriptorString +"(Number of locants disagrees with number of parents)");
				}
				for (int j = 0; j < fusionDescriptors.length; j++) {
					String fusionDescriptor = fusionDescriptors[j];
					Fragment component = multiplier>1 ? fusionComponents.get(j) : nextComponent;
					Fragment parentToUse = previousFusionLevelFragments.get(j);
					boolean simpleFusion = fusionDescriptor.split(":").length <= 1;
					if (simpleFusion){
						String[] fusionArray = determineNumericalAndLetterComponents(fusionDescriptor);
						if (fusionArray[1].length() != 0){
							int numberOfPrimes =StringTools.countTerminalPrimes(fusionArray[1]);
							if (fusionArray[0].length() == 0){
								fusionDescriptor = fusionArray[1].replaceAll("'", "");
							}
							else{
								fusionDescriptor = fusionArray[0]+ "-" +fusionArray[1].replaceAll("'", "");
							}
							if (numberOfPrimes !=j){//check the number of primes on the letter part agree with the parent to use e.g.[4,5-bcd:1,2-c']difuran
								throw new StructureBuildingException("Incorrect number of primes in fusion descriptor: " + fusionDescriptor);
							}
						}
						performSimpleFusion(fusionDescriptor, component, parentToUse);
					}
					else{
						// MolLangData: we should have a debug warning here, we should take care when this is happening
						state.addWarning(OpsinWarning.OpsinWarningType.MolLangData_DEBUG_WARNING, "MolLangData does not know when the code can run into this situation for this higher order fusion; please report this issue to the developers");
						performHigherOrderFusion(fusionDescriptor, component, parentToUse);
					}
					// MolLangData TODO: we should merge the same component into a single
				}
				previousFusionLevelFragments = fusionComponents;
				componentFragments.addAll(fusionComponents);
			}
			if (previousFusionLevelFragments.size()!=1){
				throw new StructureBuildingException("Invalid fused ring system. Incomplete multiparent system");
			}
		}
		return i;
	}

	/**
	 * Splits a first order fusion component into it's numerical and letter parts
	 * Either one of these can be the blank string as they may have been omitted
	 * The first entry in the array is the numbers and the second the letters
	 * @param fusionDescriptor
	 * @return
	 */
	private String[] determineNumericalAndLetterComponents(String fusionDescriptor) {
		String[] fusionArray = fusionDescriptor.split("-");
		if (fusionArray.length ==2){
			return fusionArray;
		}
		else{
			String[] components = new String[2];
			if (fusionArray[0].contains(",")){//the digit section
				components[0]=fusionArray[0];
				components[1]="";
			}
			else{
				components[0]="";
				components[1]=fusionArray[0];
			}
			return components;
		}
	}

	/**
	 * Searches groups for benz(o) components and fuses them in accordance with
	 * FR-2.2.8 Heterobicyclic components with a benzene ring
	 * For outputParse, we do NOT detach or remove groups
	 * @throws StructureBuildingException
	 */
	private void processBenzoFusions() throws StructureBuildingException {
		for(int i = groupsInFusedRing.size() - 2; i >= 0; i--) {
			Element group = groupsInFusedRing.get(i);
			if (group.getValue().equals("benz") || group.getValue().equals("benzo")) {
				Element possibleFusionbracket = OpsinTools.getNextSibling(group);
				if (!possibleFusionbracket.getName().equals(FUSION_EL)) {
					Element possibleMultiplier = OpsinTools.getPreviousSibling(group);
					if (possibleMultiplier == null || !possibleMultiplier.getName().equals(MULTIPLIER_EL) || possibleMultiplier.getAttributeValue(TYPE_ATR).equals(GROUP_TYPE_VAL)) {
						//e.g. 2-benzofuran. Fused rings of this type are a special case treated as being a single component
						//and have a special convention for indicating the position of heteroatoms 
						benzoSpecificFusion(group, groupsInFusedRing.get(i + 1));
						group.detach();
						groupsInFusedRing.remove(i);
					}
				}
			}
		}
	}

	/**
	 * Modifies nextComponent's locants according to the fusionLevel.
	 * @param component
	 * @param fusionLevel
	 */
	private void relabelAccordingToFusionLevel(Fragment component, int fusionLevel)  {
		if (fusionLevel > 0){
			FragmentTools.relabelNumericLocants(component.getAtomList(), StringTools.multiplyString("'", fusionLevel));
		}
	}

	/**
	 * Handles fusion between components where the fusion descriptor is of the form:
	 * comma separated locants dash letters
	 * e.g imidazo[4,5-d]pyridine
	 * The fusionDescriptor may be given as null or the letter/numerical part omitted.
	 * Sensible defaults will be found instead
	 * @param fusionDescriptor
	 * @param childRing
	 * @param parentRing
	 * @throws StructureBuildingException
	 */
	private void performSimpleFusion(String fusionDescriptor, Fragment childRing, Fragment parentRing) throws StructureBuildingException {
		List<String> numericalLocantsOfChild = null;
		List<String> letterLocantsOfParent = null;
		if (fusionDescriptor != null){
			String[] fusionArray = fusionDescriptor.split("-");
			if (fusionArray.length ==2){
				numericalLocantsOfChild = Arrays.asList(fusionArray[0].split(","));
				char[] tempLetterLocantsOfParent = fusionArray[1].toCharArray();
				letterLocantsOfParent = new ArrayList<>();
                for (char letterLocantOfParent : tempLetterLocantsOfParent) {
                    letterLocantsOfParent.add(String.valueOf(letterLocantOfParent));
                }
			}
			else{
				if (fusionArray[0].contains(",")){//only has digits
					String[] numericalLocantsOfChildTemp = fusionArray[0].split(",");
					numericalLocantsOfChild = Arrays.asList(numericalLocantsOfChildTemp);
				}
				else{//only has letters
					char[] tempLetterLocantsOfParentCharArray = fusionArray[0].toCharArray();
					letterLocantsOfParent = new ArrayList<>();
                    for (char letterLocantOfParentCharArray : tempLetterLocantsOfParentCharArray) {
                        letterLocantsOfParent.add(String.valueOf(letterLocantOfParentCharArray));
                    }
				}
			}
		}

		int edgeLength =1;
		if (numericalLocantsOfChild != null){
			if (numericalLocantsOfChild.size() <=1){
				throw new StructureBuildingException("At least two numerical locants must be provided to perform fusion!");
			}
			edgeLength = numericalLocantsOfChild.size()-1;
		}
		else if (letterLocantsOfParent != null){
			edgeLength = letterLocantsOfParent.size();
		}

		if (numericalLocantsOfChild == null){
			numericalLocantsOfChild = findPossibleNumericalLocants(childRing, edgeLength);
		}

		if (letterLocantsOfParent == null){
			letterLocantsOfParent = findPossibleLetterLocants(parentRing, edgeLength);
		}
		if (numericalLocantsOfChild == null || letterLocantsOfParent ==null){
			throw new StructureBuildingException("Unable to find bond to form fused ring system. Some information for forming fused ring system was only supplyed implicitly");
		}

		processFirstOrderFusionDescriptors(childRing, parentRing, numericalLocantsOfChild, letterLocantsOfParent);//fuse the rings
	}

	/**
	 * Takes a ring an returns and array with one letter corresponding to a side/s
	 * that contains two adjacent non bridgehead carbons
	 * The number of sides is specified by edgeLength
	 * @param ring
	 * @param edgeLength The number of bonds to be fused along
	 * @return
	 */
	private List<String> findPossibleLetterLocants(Fragment ring, int edgeLength) {
		List<Integer> carbonAtomIndexes = new ArrayList<>();
		int numberOfAtoms = ring.getAtomCount();
		CyclicAtomList cyclicAtomList = new CyclicAtomList(ring.getAtomList());
		for (int i = 0; i <= numberOfAtoms; i++) {
			//iterate backwards in list to use highest locanted edge in preference.
			//this retains what is currently locant 1 on the parent ring as locant 1 if the first two atoms found match
			//the last atom in the list is potentially tested twice e.g. on a 6 membered ring, 6-5 and 1-6 are both possible
			Atom atom = cyclicAtomList.previous();
			//want non-bridgehead carbon atoms. Double-check that these carbon atoms are actually bonded (e.g. von baeyer systems have non-consecutive atom numbering!)
			if (atom.getElement() == ChemEl.C && atom.getBondCount() == 2
					&& (carbonAtomIndexes.isEmpty() || atom.getAtomNeighbours().contains(cyclicAtomList.peekNext()))){
				carbonAtomIndexes.add(cyclicAtomList.getIndex());
				if (carbonAtomIndexes.size() == edgeLength + 1){//as many carbons in a row as to give that edgelength ->use these side/s
					Collections.reverse(carbonAtomIndexes);
					List<String> letterLocantsOfParent = new ArrayList<>();
					for (int j = 0; j < edgeLength; j++) {
						letterLocantsOfParent.add(String.valueOf((char)(97 + carbonAtomIndexes.get(j))));//97 is ascii for a	
					}
					return letterLocantsOfParent;
				}
			}
			else{
				carbonAtomIndexes.clear();
			}
		}
		return null;
	}

	/**
	 * Takes a ring and returns an array of numbers corresponding to a side/s
	 * that contains two adjacent non bridgehead carbons
	 * The number of sides is specified by edgeLength
	 * @param ring
	 * @param edgeLength The number of bonds to be fused along
	 * @return
	 */
	private List<String> findPossibleNumericalLocants(Fragment ring, int edgeLength) {
		List<String> carbonLocants = new ArrayList<>();
		int numberOfAtoms = ring.getAtomCount();
		CyclicAtomList cyclicAtomList = new CyclicAtomList(ring.getAtomList());
		for (int i = 0; i <= numberOfAtoms; i++) {
			//the last atom in the list is potentially tested twice e.g. on a 6 membered ring, 1-2 and 6-1 are both possible
			Atom atom = cyclicAtomList.next();
			//want non-bridgehead carbon atoms. Double-check that these carbon atoms are actually bonded (e.g. von baeyer systems have non-consecutive atom numbering!)
			if (atom.getElement() == ChemEl.C && atom.getBondCount() == 2
					&& (carbonLocants.isEmpty() || atom.getAtomNeighbours().contains(cyclicAtomList.peekPrevious()))){
				carbonLocants.add(atom.getFirstLocant());
				if (carbonLocants.size() == edgeLength + 1){//as many carbons in a row as to give that edgelength ->use these side/s
					List<String> numericalLocantsOfChild = new ArrayList<>();
					for (String locant : carbonLocants) {
						numericalLocantsOfChild.add(locant);
					}
					return numericalLocantsOfChild;
				}
			}
			else{
				carbonLocants.clear();
			}
		}
		return null;
	}

	/**
	 * Performs a single ring fusion using the values in numericalLocantsOfChild/letterLocantsOfParent
	 * @param childRing
	 * @param parentRing
	 * @param numericalLocantsOfChild
	 * @param letterLocantsOfParent
	 * @throws StructureBuildingException
	 */
	private void processFirstOrderFusionDescriptors(Fragment childRing, Fragment parentRing, List<String> numericalLocantsOfChild, List<String> letterLocantsOfParent) throws StructureBuildingException {
		List<Atom> childAtoms = determineAtomsToFuse(childRing, numericalLocantsOfChild, letterLocantsOfParent.size() +1);
		if (childAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}

		List<Atom> parentAtoms = new ArrayList<>();
		List<Atom> parentPeripheralAtomList = getPeripheralAtoms(parentRing.getAtomList());
		CyclicAtomList cyclicListAtomsOnSurfaceOfParent = new CyclicAtomList(parentPeripheralAtomList, (int)letterLocantsOfParent.get(0).charAt(0) -97);//convert from lower case character through ascii to 0-23
		parentAtoms.add(cyclicListAtomsOnSurfaceOfParent.getCurrent());
		for (int i = 0; i < letterLocantsOfParent.size(); i++) {
			parentAtoms.add(cyclicListAtomsOnSurfaceOfParent.next());
		}
		fuseRings(childAtoms, parentAtoms);
		// MolLangData: merge the child ring tokenEl to the parent ring tokenEl by creating a FusedChildRing element
		// mergeChildRingTokenElToParentRingTokenEl(childRing, parentRing, childAtoms, parentAtoms);
	}
	
	/**
	 * Merges the childRing tokenEl to the parentRing tokenEl by creating a FusedChildRing element
	 * that contains fusion information about which atoms in parent and child are fused.
	 * @param childRing The child ring fragment
	 * @param parentRing The parent ring fragment
	 * @param childAtoms The atoms in the child ring that are fused
	 * @param parentAtoms The atoms in the parent ring that are fused
	 * @throws StructureBuildingException
	 */
	private void mergeChildRingTokenElToParentRingTokenEl(Fragment childRing, Fragment parentRing, List<Atom> childAtoms, List<Atom> parentAtoms) throws StructureBuildingException {
		// MolLangData: directly raise an exception here, because this is deprecated and should not be used
		throw new UnsupportedOperationException("mergeChildRingTokenElToParentRingTokenEl is deprecated and should not be used");
		
		/* 
		Element childRingTokenEl = childRing.getTokenEl();
		Element parentRingTokenEl = parentRing.getTokenEl();
		
		if (childRingTokenEl == null || parentRingTokenEl == null) {
			return; // Cannot merge if tokenEls are not available
		}
		// create a new FusedChildRing to store the 


		
		// Get atom lists in the order they appear in atomMapFromId (which corresponds to labels order)
		List<Atom> childAtomList = childRing.getAtomList();
		List<Atom> parentAtomList = parentRing.getAtomList();
		
		// Get labels from tokenEls
		String childLabelsStr = childRingTokenEl.getAttributeValue(LABELS_ATR);
		String parentLabelsStr = parentRingTokenEl.getAttributeValue(LABELS_ATR);
		
		if (childLabelsStr == null || parentLabelsStr == null) {
			return; // Cannot determine fusion labels if labels are not available
		}
		
		String[] childLabels = childLabelsStr.split("/", -1);
		String[] parentLabels = parentLabelsStr.split("/", -1);
		
		// Find the label indices for fused atoms
		List<String> fusedChildLabels = new ArrayList<>();
		List<String> fusedParentLabels = new ArrayList<>();
		
		for (Atom childAtom : childAtoms) {
			int index = childAtomList.indexOf(childAtom);
			if (index >= 0 && index < childLabels.length) {
				fusedChildLabels.add(childLabels[index]);
			}
		}
		
		for (Atom parentAtom : parentAtoms) {
			int index = parentAtomList.indexOf(parentAtom);
			if (index >= 0 && index < parentLabels.length) {
				fusedParentLabels.add(parentLabels[index]);
			}
		}
		
		// Create FusedChildRing element
		Element fusedChildRingEl = new TokenEl(FUSEDCHILDRING_EL);
		
		// Set the value from childRing tokenEl value
		fusedChildRingEl.setValue(childRingTokenEl.getValue());
		
		// Copy all attributes from childRing tokenEl to fusedChildRingEl
		for (int i = 0; i < childRingTokenEl.getAttributeCount(); i++) {
			Attribute attr = childRingTokenEl.getAttribute(i);
			fusedChildRingEl.addAttribute(new Attribute(attr));
		}
		
		// Copy all children from childRing tokenEl to fusedChildRingEl
		for (Element childEl : childRingTokenEl.getChildElements()) {
			Element copiedChild = childEl.copy();
			copiedChild.setParent(fusedChildRingEl);
			fusedChildRingEl.addChild(copiedChild);
		}
		
		//fusedChildRingEl.addAttribute(new Attribute("fusedChildLabels", fusedChildLabelsStr.toString()));
		//fusedChildRingEl.addAttribute(new Attribute("fusedParentLabels", fusedParentLabelsStr.toString()));
		
		// Insert FusedChildRing element after parentRing tokenEl
		// Directly combine the parentRing tokenEl and the fusedChildRingEl
		parentRingTokenEl.addChild(fusedChildRingEl);

		*/
	}
	
	/**
	 * Returns the sublist of the given atoms that are peripheral atoms given that the list is ordered such that the interior atoms are at the end of the list
	 * @param atomList
	 * @return
	 */
	private List<Atom> getPeripheralAtoms(List<Atom> atomList) {
		//find the indice of the last atom on the surface of the ring. This obviously connects to the first atom. The objective is to exclude any interior atoms.
		List<Atom> neighbours = atomList.get(0).getAtomNeighbours();
		int indice = Integer.MAX_VALUE;
		for (Atom atom : neighbours) {
			int indexOfAtom =atomList.indexOf(atom);
			if (indexOfAtom ==1){//not the next atom
				continue;
			}
			else if (indexOfAtom ==-1){//not in parentRing
				continue;
			}
			if (atomList.indexOf(atom)< indice){
				indice = indexOfAtom;
			}
		}
		return atomList.subList(0, indice +1);
	}

	/**
	 * Handles fusion between components where the fusion descriptor is of the form:
	 * comma separated locants colon comma separated locants
	 * e.g pyrido[1'',2'':1',2']imidazo
	 * @param fusionDescriptor
	 * @param nextComponent
	 * @param fusedRing
	 * @throws StructureBuildingException 
	 */
	private void performHigherOrderFusion(String fusionDescriptor, Fragment nextComponent, Fragment fusedRing) throws StructureBuildingException {
		List<String> numericalLocantsOfChild = null;
		List<String> numericalLocantsOfParent = null;
		String[] fusionArray = fusionDescriptor.split(":");
		if (fusionArray.length ==2){
			numericalLocantsOfChild = Arrays.asList(fusionArray[0].split(","));
			numericalLocantsOfParent = Arrays.asList(fusionArray[1].split(","));
		}
		else{
			throw new StructureBuildingException("Malformed fusion bracket: This is an OPSIN bug, check regexTokens.xml");
		}
		processHigherOrderFusionDescriptors(nextComponent, fusedRing, numericalLocantsOfChild, numericalLocantsOfParent);//fuse the rings
	}

	/**
	 * Performs a single ring fusion using the values in numericalLocantsOfChild/numericalLocantsOfParent
	 * @param childRing
	 * @param parentRing
	 * @param numericalLocantsOfChild
	 * @param numericalLocantsOfParent
	 * @throws StructureBuildingException
	 */
	private void processHigherOrderFusionDescriptors(Fragment childRing, Fragment parentRing, List<String> numericalLocantsOfChild, List<String> numericalLocantsOfParent) throws StructureBuildingException {
		List<Atom> childAtoms =determineAtomsToFuse(childRing, numericalLocantsOfChild, null);
		if (childAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}

		List<Atom> parentAtoms = determineAtomsToFuse(parentRing, numericalLocantsOfParent, childAtoms.size());
		if (parentAtoms ==null){
			throw new StructureBuildingException("Malformed fusion bracket!");
		}
		fuseRings(childAtoms, parentAtoms);
	}

	/**
	 * Determines which atoms on a ring should be used for fusion given a set of numerical locants.
	 * If from the other ring involved in the fusion it is known how many atoms are expected to be found this should be provided
	 * If this is not known it should be set to null and the smallest number of fusion atoms will be returned.
	 * @param ring
	 * @param numericalLocantsOnRing
	 * @param expectedNumberOfAtomsToBeUsedForFusion
	 * @return
	 * @throws StructureBuildingException
	 */
	private List<Atom> determineAtomsToFuse(Fragment ring, List<String> numericalLocantsOnRing, Integer expectedNumberOfAtomsToBeUsedForFusion) throws StructureBuildingException {
		List<Atom> parentPeripheralAtomList = getPeripheralAtoms(ring.getAtomList());
		String firstLocant = numericalLocantsOnRing.get(0);
		String lastLocant = numericalLocantsOnRing.get(numericalLocantsOnRing.size() - 1);
		int indexfirst = parentPeripheralAtomList.indexOf(ring.getAtomByLocantOrThrow(firstLocant));
		if (indexfirst == -1) {
			throw new StructureBuildingException(firstLocant + " refers to an atom that is not a peripheral atom!");
		}
		int indexfinal = parentPeripheralAtomList.indexOf(ring.getAtomByLocantOrThrow(lastLocant));
		if (indexfinal == -1) {
			throw new StructureBuildingException(lastLocant + " refers to an atom that is not a peripheral atom!");
		}
		CyclicAtomList cyclicRingAtomList = new CyclicAtomList(parentPeripheralAtomList, indexfirst);
		List<Atom> fusionAtoms = null;
		
		List<Atom> potentialFusionAtomsAscending = new ArrayList<>();
		potentialFusionAtomsAscending.add(cyclicRingAtomList.getCurrent());
		while (cyclicRingAtomList.getIndex() != indexfinal){//assume numbers are ascending
			potentialFusionAtomsAscending.add(cyclicRingAtomList.next());
		}
		if (expectedNumberOfAtomsToBeUsedForFusion ==null ||expectedNumberOfAtomsToBeUsedForFusion == potentialFusionAtomsAscending.size()){
			boolean notInPotentialParentAtoms =false;
			for (int i =1; i < numericalLocantsOnRing.size()-1 ; i ++){
				if (!potentialFusionAtomsAscending.contains(ring.getAtomByLocantOrThrow(numericalLocantsOnRing.get(i)))){
					notInPotentialParentAtoms =true;
				}
			}
			if (!notInPotentialParentAtoms){
				fusionAtoms = potentialFusionAtomsAscending;
			}
		}
		
		if (fusionAtoms ==null || expectedNumberOfAtomsToBeUsedForFusion ==null){//that didn't work, so try assuming the numbers are descending
			cyclicRingAtomList.setIndex(indexfirst);
			List<Atom> potentialFusionAtomsDescending = new ArrayList<>();
			potentialFusionAtomsDescending.add(cyclicRingAtomList.getCurrent());
			while (cyclicRingAtomList.getIndex() != indexfinal){//assume numbers are descending
				potentialFusionAtomsDescending.add(cyclicRingAtomList.previous());
			}
			if (expectedNumberOfAtomsToBeUsedForFusion ==null || expectedNumberOfAtomsToBeUsedForFusion == potentialFusionAtomsDescending.size()){
				boolean notInPotentialParentAtoms =false;
				for (int i =1; i < numericalLocantsOnRing.size()-1 ; i ++){
					if (!potentialFusionAtomsDescending.contains(ring.getAtomByLocantOrThrow(numericalLocantsOnRing.get(i)))){
						notInPotentialParentAtoms =true;
					}
				}
				if (!notInPotentialParentAtoms){
					if (fusionAtoms!=null && expectedNumberOfAtomsToBeUsedForFusion ==null){
						//prefer less fusion atoms
						if (potentialFusionAtomsDescending.size()< fusionAtoms.size()){
							fusionAtoms = potentialFusionAtomsDescending;
						}
					}
					else{
						fusionAtoms = potentialFusionAtomsDescending;
					}
				}
			}
		}
		return fusionAtoms;
	}

	/**
	 * Creates the bonds required to fuse two rings together.
	 * The child atoms are recorded as atoms that should be removed later
	 * @param childAtoms
	 * @param parentAtoms
	 * @throws StructureBuildingException
	 */
	private void fuseRings(List<Atom> childAtoms, List<Atom> parentAtoms) throws StructureBuildingException {
		if (parentAtoms.size()!=childAtoms.size()){
			throw new StructureBuildingException("Problem with fusion descriptors: Parent atoms specified: " + parentAtoms.size() +" Child atoms specified: " + childAtoms.size() + " These should have been identical!");
		}
		//replace parent atoms if the atom has already been used in fusion with the original atom
		//This will occur if fusion has resulted in something resembling a spiro centre e.g. cyclopenta[1,2-b:5,1-b']bis[1,4]oxathiine
		for (int i = parentAtoms.size() -1; i >=0; i--) {
			if (atomsToRemoveToReplacementAtom.get(parentAtoms.get(i))!=null){
				parentAtoms.set(i, atomsToRemoveToReplacementAtom.get(parentAtoms.get(i)));
			}
			if (atomsToRemoveToReplacementAtom.get(childAtoms.get(i))!=null){
				childAtoms.set(i, atomsToRemoveToReplacementAtom.get(childAtoms.get(i)));
			}
		}
		
		//sync spareValency and check that element type matches
		for (int i = 0; i < childAtoms.size(); i++) {
			Atom parentAtom = parentAtoms.get(i);
			Atom childAtom = childAtoms.get(i);
			if (childAtom.hasSpareValency()){
				parentAtom.setSpareValency(true);
			}
			if (parentAtom.getElement() != childAtom.getElement()){
				throw new StructureBuildingException("Invalid fusion descriptor: Heteroatom placement is ambiguous as it is not present in both components of the fusion");
			}
			atomsToRemoveToReplacementAtom.put(childAtom, parentAtom);
		}
		
		Set<Bond> fusionEdgeBonds  = new HashSet<>();//these bonds already exist in both the child and parent atoms
		for (int i = 0; i < childAtoms.size() -1; i++) {
			fusionEdgeBonds.add(childAtoms.get(i).getBondToAtomOrThrow(childAtoms.get(i+1)));
			fusionEdgeBonds.add(parentAtoms.get(i).getBondToAtomOrThrow(parentAtoms.get(i+1)));
		}

		Set<Bond> bondsToAddToParentAtoms = new LinkedHashSet<>();
		for (Atom childAtom : childAtoms) {
			for (Bond b : childAtom.getBonds()) {
				if (!fusionEdgeBonds.contains(b)){
					bondsToAddToParentAtoms.add(b);
				}
			}
		}
		
		Set<Bond> bondsToAddToChildAtoms = new LinkedHashSet<>();
		for (Atom parentAtom : parentAtoms) {
			for (Bond b : parentAtom.getBonds()) {
				if (!fusionEdgeBonds.contains(b)){
					bondsToAddToChildAtoms.add(b);
				}
			}
		}
		
		for (Bond bond : bondsToAddToParentAtoms) {
			Atom from = bond.getFromAtom();
			int indiceInChildAtoms = childAtoms.indexOf(from);
			if (indiceInChildAtoms !=-1){
				from = parentAtoms.get(indiceInChildAtoms);
			}
			Atom to = bond.getToAtom();
			indiceInChildAtoms = childAtoms.indexOf(to);
			if (indiceInChildAtoms !=-1){
				to = parentAtoms.get(indiceInChildAtoms);
			}
			state.fragManager.createBond(from, to, 1);
		}

		for (Bond bond : bondsToAddToChildAtoms) {
			Atom from = bond.getFromAtom();
			int indiceInParentAtoms = parentAtoms.indexOf(from);
			if (indiceInParentAtoms !=-1){
				from = childAtoms.get(indiceInParentAtoms);
			}
			Atom to = bond.getToAtom();
			indiceInParentAtoms = parentAtoms.indexOf(to);
			if (indiceInParentAtoms !=-1){
				to = childAtoms.get(indiceInParentAtoms);
			}
			Bond newBond = new Bond(from, to, 1);
			if (childAtoms.contains(from)){
				from.addBond(newBond);
			}
			else{
				to.addBond(newBond);
			}
		}
	}

	/**
	 * Fuse the benzo with the subsequent ring
	 * Uses locants in front of the benz/benzo group to assign heteroatoms on the now numbered fused ring system
	 * @param benzoEl
	 * @param parentEl
	 * @throws StructureBuildingException
	 */
	private void benzoSpecificFusion(Element benzoEl, Element parentEl) throws StructureBuildingException {
		/*
		 * Perform the fusion, number it and associate it with the parentEl
		 */
		Fragment benzoRing = benzoEl.getFrag();
		Fragment parentRing = parentEl.getFrag();
		
		// MolLangData: Create snapshots of labels and atoms before fusion
		Map<Atom, String> parentRingAtomToLabel = createAtomToLabelSnapshot(parentRing);
		Map<Atom, String> fusedRingAtomToLabel = createAtomToLabelSnapshot(benzoRing);
		List<Map<Atom, String>> fusedRingAtomToLabelList = new ArrayList<>();
		fusedRingAtomToLabelList.add(fusedRingAtomToLabel);
		
		performSimpleFusion(null, benzoRing , parentRing);
		state.fragManager.incorporateFragment(benzoRing, parentRing);
		
		// MolLangData: Update snapshots to map removed atoms to their replacement atoms
		updateSnapshotWithAtomReplacements(parentRingAtomToLabel, atomsToRemoveToReplacementAtom);
		updateSnapshotWithAtomReplacements(fusedRingAtomToLabel, atomsToRemoveToReplacementAtom);
		
		removeMergedAtoms();
		FusedRingNumberer.numberFusedRing(parentRing);//numbers the fused ring;
		// MolLangData: create a fusedRingNumbering element to store the numbering information
		createFusedRingNumberingElement(parentRing, parentRingAtomToLabel, fusedRingAtomToLabelList);
		Fragment fusedRing =parentRing;
		setBenzoHeteroatomPositioning(benzoEl, fusedRing);
	}

	/**
	 * Checks for locant(s) before benzo and uses these to set 
	 * @param benzoEl
	 * @param fusedRing
	 * @throws StructureBuildingException
	 */
	private void setBenzoHeteroatomPositioning(Element benzoEl, Fragment fusedRing) throws StructureBuildingException {
		Element locantEl = OpsinTools.getPreviousSibling(benzoEl);
		if (locantEl != null && locantEl.getName().equals(LOCANT_EL)) {
			String[] locants = locantEl.getValue().split(",");
			if (locantsCouldApplyToHeteroatomPositions(locants, benzoEl)) {
				List<Atom> atomList =fusedRing.getAtomList();
				List<Atom> heteroatoms = new ArrayList<>();
				List<ChemEl> elementOfHeteroAtom = new ArrayList<>();
				for (Atom atom : atomList) {//this iterates in the same order as the numbering system
					if (atom.getElement() != ChemEl.C){
						heteroatoms.add(atom);
						elementOfHeteroAtom.add(atom.getElement());
					}
				}
				if (locants.length == heteroatoms.size()){//as many locants as there are heteroatoms to assign
					//check for special case of a single locant indicating where the group substitutes e.g. 4-benzofuran-2-yl
					if (!(locants.length == 1 && OpsinTools.getPreviousSibling(locantEl) == null
							 && ComponentProcessor.checkLocantPresentOnPotentialRoot(state, benzoEl.getParent(), locants[0]))) {
						for (Atom atom : heteroatoms) {
							atom.setElement(ChemEl.C);
						}
						for (int i=0; i< heteroatoms.size(); i++) {
							fusedRing.getAtomByLocantOrThrow(locants[i]).setElement(elementOfHeteroAtom.get(i));
						}
						locantEl.detach();
					}
				}
				else if (locants.length > 1){
					throw new StructureBuildingException("Unable to assign all locants to benzo-fused ring or multiplier was mising");
				}
			}
		}
	}

	private boolean locantsCouldApplyToHeteroatomPositions(String[] locants, Element benzoEl) {
		if (!locantsAreAllNumeric(locants)) {
			return false;
		}
		List<Element> suffixes = benzoEl.getParent().getChildElements(SUFFIX_EL);
		int suffixesWithoutLocants = 0;
		for (Element suffix : suffixes) {
			if (suffix.getAttribute(LOCANT_ATR)==null){
				suffixesWithoutLocants++;
			}
		}
		if (locants.length == suffixesWithoutLocants){//In preference locants will be assigned to suffixes rather than to this nomenclature
			return false;
		}
		return true;
	}

	private boolean locantsAreAllNumeric(String[] locants) {
		for (String locant : locants) {
			if (!MATCH_NUMERIC_LOCANT.matcher(locant).matches()){
				return false;
			}
		}
		return true;
	}

	/**
	 * MolLangData: Creates a snapshot mapping atoms to their labels
	 * First tries to get labels from atomMapFromLocant, then falls back to index-based approach from labels attribute
	 * @param ring The fragment to create snapshot for
	 * @return Map from Atom to label string
	 */
	private Map<Atom, String> createAtomToLabelSnapshot(Fragment ring) {
		Map<Atom, String> atomToLabel = new HashMap<>();
		Element tokenEl = ring.getTokenEl();
		
		if (tokenEl == null) {
			return atomToLabel;
		}
		
		List<Atom> atomList = ring.getAtomList();
		
		// First, try to get labels from atomMapFromLocant
		Set<String> locants = ring.getLocants();
		for (String locant : locants) {
			Atom atom = ring.getAtomByLocant(locant);
			if (atom != null && atomList.contains(atom)) {
				// Get the first locant from the label (in case of multiple comma-separated locants)
				String firstLabel = locant.split(",")[0].trim();
				if (!firstLabel.isEmpty()) {
					// Only set if not already set (prefer first locant found)
					if (!atomToLabel.containsKey(atom)) {
						atomToLabel.put(atom, firstLabel);
					}
				}
			}
		}
		
		// Fall back to index-based approach from labels attribute for atoms not found in atomMapFromLocant
		String labelsStr = tokenEl.getAttributeValue(LABELS_ATR);
		if (labelsStr != null) {
			String[] labels = labelsStr.split("/", -1);
			for (int i = 0; i < Math.min(labels.length, atomList.size()); i++) {
				Atom atom = atomList.get(i);
				if (!atomToLabel.containsKey(atom)) {
					String label = labels[i];
					if (label != null && !label.isEmpty()) {
						// Get the first locant from the label (in case of multiple comma-separated locants)
						String firstLabel = label.split(",")[0].trim();
						if (!firstLabel.isEmpty()) {
							atomToLabel.put(atom, firstLabel);
						}
					}
				}
			}
		}
		
		return atomToLabel;
	}

	/**
	 * MolLangData: Updates a snapshot map to replace removed atoms with their replacement atoms
	 * Replaces the key (atom) in the snapshot map with the replacement atom from atomReplacements
	 * @param snapshot The snapshot map to update
	 * @param atomReplacements The map of removed atoms to their replacement atoms
	 */
	private void updateSnapshotWithAtomReplacements(Map<Atom, String> snapshot, Map<Atom, Atom> atomReplacements) {
		// Create a list of entries to update to avoid concurrent modification
		List<Map.Entry<Atom, String>> entriesToUpdate = new ArrayList<>();
		for (Map.Entry<Atom, String> entry : snapshot.entrySet()) {
			Atom atom = entry.getKey();
			if (atomReplacements.containsKey(atom)) {
				entriesToUpdate.add(entry);
			}
		}
		
		// Update entries: replace the key (atom) with the replacement atom
		for (Map.Entry<Atom, String> entry : entriesToUpdate) {
			Atom oldAtom = entry.getKey();
			String label = entry.getValue();
			Atom replacementAtom = atomReplacements.get(oldAtom);
			// Follow the replacement chain to find the final atom
			Atom finalAtom = replacementAtom;
			while (atomReplacements.containsKey(finalAtom)) {
				finalAtom = atomReplacements.get(finalAtom);
			}
			// Remove old entry and add with replacement atom as key
			snapshot.remove(oldAtom);
			// If the final atom already exists, merge labels (comma-separated)
			if (snapshot.containsKey(finalAtom)) {
				String existingLabel = snapshot.get(finalAtom);
				snapshot.put(finalAtom, existingLabel + "," + label);
			} else {
				snapshot.put(finalAtom, label);
			}
		}
	}

	/**
	 * MolLangData: Updates atomToOriginalLabels to map removed atoms to their replacement atoms
	 * @param atomReplacements The map of removed atoms to their replacement atoms
	 */
	private void updateAtomToOriginalLabelsWithReplacements(Map<Atom, Atom> atomReplacements) {
		Map<Atom, Map<Integer, String>> updatedMap = new HashMap<>();
		for (Map.Entry<Atom, Map<Integer, String>> entry : atomToOriginalLabels.entrySet()) {
			Atom atom = entry.getKey();
			Map<Integer, String> ringLabels = entry.getValue();
			// Follow the replacement chain to find the final atom
			Atom finalAtom = atom;
			while (atomReplacements.containsKey(finalAtom)) {
				finalAtom = atomReplacements.get(finalAtom);
			}
			// Merge labels if the final atom already exists
			if (updatedMap.containsKey(finalAtom)) {
				updatedMap.get(finalAtom).putAll(ringLabels);
			} else {
				updatedMap.put(finalAtom, new HashMap<>(ringLabels));
			}
		}
		atomToOriginalLabels.clear();
		atomToOriginalLabels.putAll(updatedMap);
	}

	/**
	 * MolLangData: Creates a fusedRingNumbering element to store the numbering information
	 * of the fused ring system. The element contains:
	 * - labels: The new locants from atomMapFromLocant in parentRing (e.g., "1/2/3/3a/4/5/6/7/7a")
	 * - originalLabels: For each new label, the original label from each ring (e.g., "(1, )/(2, )/(3, )/(4,1,2)/(, 6,)/(, 5,)/(, 4,)/(,3,)/(5,2,)")
	 *   Format: (parent, fusedRing1, fusedRing2, ...)
	 * 
	 * @param fusedRing The fused ring fragment after numbering
	 * @param parentRingAtomToLabel Snapshot of parent ring atom-to-label mapping (null for general fusion)
	 * @param fusedRingAtomToLabelList List of snapshots of rings being fused into the parent ring atom-to-label mappings (null for general fusion)
	 */
	private void createFusedRingNumberingElement(Fragment fusedRing, Map<Atom, String> parentRingAtomToLabel, List<Map<Atom, String>> fusedRingAtomToLabelList) {
		Element parentRingTokenEl = fusedRing.getTokenEl();
		if (parentRingTokenEl == null) {
			return;
		}
		
		List<Atom> atomList = fusedRing.getAtomList();
		if (atomList.isEmpty()) {
			return;
		}
		
		// If fusedRingAtomToLabelList is empty, return early
		if (fusedRingAtomToLabelList != null && fusedRingAtomToLabelList.isEmpty()) {
			return;
		}
		
		// Build labels string from the new locants (ordered by atomList)
		StringBuilder labelsBuilder = new StringBuilder();
		StringBuilder originalLabelsBuilder = new StringBuilder();
		
		// Iterate through atoms in the fused ring
		for (int i = 0; i < atomList.size(); i++) {
			Atom atom = atomList.get(i);
			
			// Get the new label (first locant of the atom)
			String newLabel = atom.getFirstLocant();
			if (newLabel == null || newLabel.isEmpty()) {
				newLabel = "";
			}
			
			if (i > 0) {
				labelsBuilder.append("/");
				originalLabelsBuilder.append("/");
			}
			labelsBuilder.append(newLabel);
			
			// Get original labels from snapshots (for specific fusion) or atomToOriginalLabels (for general fusion)
			String parentLabel = "";
			List<String> childLabels = new ArrayList<>();
			
			if (parentRingAtomToLabel != null && fusedRingAtomToLabelList != null && !fusedRingAtomToLabelList.isEmpty()) {
				// Specific fusion case with snapshots: use snapshots
				parentLabel = parentRingAtomToLabel.getOrDefault(atom, "");
				for (Map<Atom, String> fusedRingAtomToLabel : fusedRingAtomToLabelList) {
					String label = fusedRingAtomToLabel.getOrDefault(atom, "");
					childLabels.add(label);
				}
			} else {
				// General fusion case: use atomToOriginalLabels
				Map<Integer, String> ringLabels = atomToOriginalLabels.get(atom);
				if (ringLabels != null) {
					int parentRingIndex = groupsInFusedRing.size() - 1;
					for (Map.Entry<Integer, String> entry : ringLabels.entrySet()) {
						Integer ringIndex = entry.getKey();
						String label = entry.getValue();
						if (label != null && !label.isEmpty()) {
							if (ringIndex == parentRingIndex) {
								parentLabel = label;
							} else {
								// For child rings, collect all labels
								childLabels.add(label);
							}
						}
					}
				}
				
				// If still not found, try to get from parent ring's tokenEl
				if (parentLabel.isEmpty()) {
					int parentRingIndex = groupsInFusedRing.size() - 1;
					Element parentGroup = groupsInFusedRing.get(parentRingIndex);
					Fragment parentRingFrag = parentGroup.getFrag();
					Element parentTokenEl = parentRingFrag.getTokenEl();
					if (parentTokenEl != null) {
						String parentLabelsStr = parentTokenEl.getAttributeValue(LABELS_ATR);
						if (parentLabelsStr != null) {
							String[] parentLabels = parentLabelsStr.split("/", -1);
							List<Atom> parentAtomList = parentRingFrag.getAtomList();
							int atomIndex = parentAtomList.indexOf(atom);
							if (atomIndex >= 0 && atomIndex < parentLabels.length) {
								String label = parentLabels[atomIndex].split(",")[0].trim();
								if (!label.isEmpty()) {
									parentLabel = label;
								}
							}
						}
					}
				}
			}
			
			// Format: (parent, fusedRing1, fusedRing2, ...)
			originalLabelsBuilder.append("(").append(parentLabel);
			for (String childLabel : childLabels) {
				originalLabelsBuilder.append(", ").append(childLabel);
			}
			originalLabelsBuilder.append(")");
		}
		
		// Create the fusedRingNumbering element
		Element fusedRingNumberingEl = new TokenEl("fusedRingLabels");
		fusedRingNumberingEl.addAttribute(new Attribute("labels", labelsBuilder.toString()));
		fusedRingNumberingEl.addAttribute(new Attribute("originalLabels", originalLabelsBuilder.toString()));
		
		// MolLangData: Add as child of parent ring tokenEl
		parentRingTokenEl.addChild(fusedRingNumberingEl);
	}
}

