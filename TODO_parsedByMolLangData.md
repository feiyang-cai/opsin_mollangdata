# TODO: Add parsedByMolLangData="false" to Token Elements

This is a checklist of XML resource files that need to have `parsedByMolLangData="false"` added to all `<token>` elements.

## Completed Files ✓
- [x] [alkanes.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/alkanes.xml)
- [x] [aminoAcids.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/aminoAcids.xml)
- [x] [arylGroups.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/arylGroups.xml)
- [x] [atomHydrides.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/atomHydrides.xml)
- [x] [carbohydrates.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/carbohydrates.xml)
- [x] [carbohydrateSuffixes.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/carbohydrateSuffixes.xml)
- [x] [chargeAndOxidationNumberSpecifiers.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/chargeAndOxidationNumberSpecifiers.xml)
- [x] [cyclicUnsaturableHydrocarbon.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/cyclicUnsaturableHydrocarbon.xml)
- [x] [elementaryAtoms.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/elementaryAtoms.xml)
- [x] [functionalTerms.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/functionalTerms.xml)
- [x] [heteroAtoms.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/heteroAtoms.xml)
- [x] [inlineSuffixes.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/inlineSuffixes.xml)
- [x] [germanTokens.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/germanTokens.xml)
- [x] [groupStemsAllowingAllSuffixes.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/groupStemsAllowingAllSuffixes.xml)

## Files To Check

### High Priority (User Mentioned)
- [ ] [arylSubstituents.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/arylSubstituents.xml)
- [ ] [carboxylicAcids.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/carboxylicAcids.xml)
- [ ] [fusionComponents.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/fusionComponents.xml)

### Additional Files
- [ ] [miscTokens.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/miscTokens.xml)
- [ ] [multipliers.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/multipliers.xml)
- [ ] [multiRadicalSubstituents.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/multiRadicalSubstituents.xml)
- [ ] [naturalProducts.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/naturalProducts.xml)
- [ ] [nonCarboxylicAcids.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/nonCarboxylicAcids.xml)
- [ ] [simpleCyclicGroups.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/simpleCyclicGroups.xml)
- [ ] [simpleGroups.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/simpleGroups.xml)
- [ ] [simpleSubstituents.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/simpleSubstituents.xml)
- [ ] [substituents.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/substituents.xml)
- [ ] [suffixes.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/suffixes.xml)
- [ ] [unsaturators.xml](opsin-core/src/main/resources/uk/ac/cam/ch/wwmm/opsin/resources/unsaturators.xml)

## Notes

- The attribute should be added to all `<token>` elements in the format: `parsedByMolLangData="false"`
- Do NOT add the attribute to `<tokenList>` elements or DOCTYPE declarations
- Some files may have specific tokenLists that should be excluded (like carbohydrates.xml had exclusions for carbohydrateRingSize and stereoChemistry)
- Check each file manually to ensure proper application
