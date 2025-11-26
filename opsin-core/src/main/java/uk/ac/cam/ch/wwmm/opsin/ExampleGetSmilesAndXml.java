package uk.ac.cam.ch.wwmm.opsin;

/**
 * Example class to demonstrate getting SMILES and XML for a chemical name.
 * This example parses a chemical name and shows how to access both SMILES and XML.
 * 
 * To run this example:
 *   javac -cp "target/classes:target/dependency/*" ExampleGetSmilesAndXml.java
 *   java -cp ".:target/classes:target/dependency/*" uk.ac.cam.ch.wwmm.opsin.ExampleGetSmilesAndXml
 */
public class ExampleGetSmilesAndXml {
	
	public static void main(String[] args) {
		//String chemicalName = "3-[4-[(2S)-2-[(3R)-3-(fluoromethyl)pyrrolidin-1-yl]propoxy]phenoxy]-2-(4-hydroxyphenyl)-1-benzothiophen-5-ol";
		//String chemicalName = "benzothiophen-5-ol";
		//String chemicalName = "5-heptyl-5,6,7,8-tetrahydrobenzo[f]benzimidazol-2-amine";
		//String chemicalName = "5,6,7,8-tetrahydrobenzo[f]benzimidazol-2-amine";
        //String chemicalName = "benzo[e]benzimidazol-2-amine";
		//String chemicalName = "1-butyl-5-[4-(diethylamino)phenyl]imino-4-methyl-2,6-dioxopyridine-3-carbonitrile";
		String chemicalName = "(2R,3R)-2-(3,4-dihydroxyphenyl)-3,4-dihydro-2H-chromene-3,5,7-triol";

		// need to deal with the following:
		//String chemicalName = "N-[3-[2-[4-(2-methylquinolin-5-yl)piperazin-1-yl]ethyl]phenyl]pyrazine-2-carboxamide";
		System.out.println("Parsing chemical name: " + chemicalName);
		System.out.println("==========================================");
		
		// Get NameToStructure instance
		NameToStructure n2s = NameToStructure.getInstance();
		
		// Parse the chemical name and get SMILES
		OpsinResult result = n2s.parseChemicalName(chemicalName);
		
		// Display results
		System.out.println("\nStatus: " + result.getStatus());
		
		if (result.getStatus() == OpsinResult.OPSIN_RESULT_STATUS.SUCCESS) {
			String smiles = result.getSmiles();
			System.out.println("\nSMILES: " + smiles);
			
			// Get XML directly from OpsinResult
			String parseXml = result.getParseXml();
			String outputXml = result.getOutputXml();
			
			if (parseXml != null) {
				System.out.println("\nParse XML:");
				System.out.println("----------------------------------------");
				System.out.println(parseXml);
				System.out.println("----------------------------------------");
			}
			
			if (outputXml != null) {
				System.out.println("\nOutput XML:");
				System.out.println("----------------------------------------");
				System.out.println(outputXml);
				System.out.println("----------------------------------------");
			}
		} else {
			System.out.println("\nParsing failed: " + result.getMessage());
		}
	}
}

