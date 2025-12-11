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
		//String chemicalName = "(2R,3R)-2-(3,4-dihydroxyphenyl)-3,4-dihydro-2H-chromene-3,5,7-triol";
		//String chemicalName = "(3R,4R,5R)-ethyl 4-acetamido-5-amino-3-(pentan-3-yloxy)cyclohex-1-enecarboxylate";
		//String chemicalName = "(2Z)-5-ethyl-6-methylhex-2-ene";
		//String chemicalName = "cis-1,2-dimethylcyclohexane";
		//String chemicalName = "cis-1,2-dichloroethene";
		//String chemicalName = "(2E,4Z)-hexa-2,4-diene";
		//String chemicalName = "(3E)-3-methyl-N-[(3Z)-2-methylpenta-1,3-dien-3-yl]hexa-3,5-dien-2-imine";
		//String chemicalName = "spiro-9,9'-bifluoren";
		//String chemicalName = "2-phenyl-3,4-dihydro-2H-1-benzopyran";
		//String chemicalName = "1-benzopyran";
		//String chemicalName = "(6aR,11aR)-6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene";

		// different stereochemistry cases

		// TODO:
		//String chemicalName = "manno-hexopyranose";
		//String chemicalName = "alpha-D-glucopyranose";
		//String chemicalName = "beta-D-glucopyranose";

		//String chemicalName = "erythro 2,3-butanediol";
		//String chemicalName = "2-deoxy-D-erythro-pentose";
		//String chemicalName = "pentanthioylium";
		//String chemicalName = "3,4-epoxybutanol";
		//String chemicalName = "5,5'-spirobi[benzo[b]phosphindolium]";
		//String chemicalName = "(+)-D-glucose";
		String chemicalName = "(+)-2-benzyl-2-methoxybenzo[1,2-b:3,4-b′]difuran-3(2H)-one";
		
		// spiro system test
		// String chemicalName = "(2R,4S,4aS)-rel-11-fluoro-2,4-dimethyl-8-(methylsulfinyl)-1,2,4,4a-tetrahydro-2′H,6H-spiro[1,4-oxazino[4,3-a][1,2]oxazolo[4,5-g]quinoline-5,5′-pyrimidine]-2′,4′,6′(1′H,3′H)-trione";
		//String chemicalName = "spiro-9,9'-bifluoren";


		//String chemicalName = "dispiro[5.1.7.2]heptadecane";
		//String chemicalName = "pentaspiro[2.0.24.0.27.0.210.0.213.03]pentadecane";
		//String chemicalName = "spiro[3.4]octane";
		//String chemicalName = "dispiro[fluorene-9,1'-cyclohexane-4',1''-indene]";
		//String chemicalName = "1,1'-spirobiindene";
		//String chemicalName = "5lambda7,5',5''-spiroter[benzo[b]phosphindol]-5-ide";
		// This is important, please carefully check the results.
		//String chemicalName = "3,3':6',6''-dispiroter[bicyclo[3.1.0]hexane]";// bridge and spiro system
		//String chemicalName = "bicyclo[3.1.0]hexane";
		//String chemicalName = "cyclopentanespirocyclobutane";
		//String chemicalName = "hexahydroazepinium-1-spiro-1'-imidazolidine-3'-spiro-1''-piperidinium dibromide";
		//String chemicalName = "cyclopentanespirocyclobutane";
		//String chemicalName = "2-cyclohexenespiro-(2'-cyclopentene)";
		//String chemicalName = "dispiro[5.1.7.2]heptadecane";
		//String chemicalName = "spiro[3.4]octane";
		//String chemicalName = "spiro[4.5]deca-1,6-diene";
		//String chemicalName = "spiro[4.5]deca-1,6-dien-2-yl";
		//String chemicalName = "5lambda^5,5'-spirobi[benzo[b]phosphindol]-5-ylium";
		//String chemicalName = "1H-2lambda5-spiro[isoquinoline-2,2'-pyrido[1,2-a]pyrazin]-2-ylium";

		// fused / bridge system test
		//String chemicalName = "cyclopenta[1,2-b:5,1-b']bis[1,4]oxathiine";
		//String chemicalName = "1,12-ethenobenzo[4,5]cyclohepta[1,2,3-de]naphthalene";
		//String chemicalName = "12-methyl-2-azatricyclo[4.4.3.0¹,⁶]trideca-2,4,7,9-tetraene";
		//String chemicalName = "6,13-ethano-6,13-methanodibenzo[b,g][1,6]diazecine";
		//String chemicalName = "(10R)-7-amino-16-cyclopropyl-12-fluoro-2,10-dimethyl-15-oxo-10,15,16,17-tetrahydro-2H-8,4-(azeno)pyrazolo[4,3-h][2,5,11]benzoxadiazacyclotetradecine-3-carbonitrile";
		//String chemicalName = "11-chloro-9,10-(epoxymethano)anthracene";
		//String chemicalName = "2H-3,5-(epoxymethano)furo[3,4-b]pyran";


		// conjunctive nomenclature test
		//String chemicalName = "benzeneacetic acid";



		//String chemicalName = "1,4-oxazino[4,3-a][1,2]oxazolo[4,5-g]quinoline";


		//Hydro test
		//String chemicalName = "2(1H)-quinolinone";
		//String chemicalName = "(2R,3R)-2-(3,4-dihydroxyphenyl)-3,4-dihydro-2H-chromene-3,5,7-triol";
		//String chemicalName = "tetrahydrofuran";

		//unsaturator test
		//String chemicalName = "1-butene";
		//String chemicalName = "cyclohex-2-ene";
		//String chemicalName = "cyclohexene";

		//heteroatom test
		//String chemicalName = "2-azabenzofuran";
		//String chemicalName = "thiazole";

		//subtractive prefix test
		//String chemicalName = "2'-deoxyadenosine";

		//isotope specification test
		//String chemicalName = "(²H1)methane";

		//multiplicative test
		//String chemicalName = "1H,3H-imidazo[1,2-a]pyridine";

		//hantzsch-widman ring test
		//String chemicalName = "1,3-oxazol-2-one";
		//String chemicalName = "2H-1,2,3-triazole";
		//String chemicalName = "2H-oxepine";
		//String chemicalName = "oxazole";

		// additive suffix test
		//String chemicalName = "methylsulfonamidobenzene";

		// conjunctive suffix test
		//String chemicalName = "benzenemethanol";

		// aric acid test
		//String chemicalName = "2-hydroxyglutaric acid";
		// dialdose test
		//String chemicalName = "L-threo-Tetrodialdose";
		// diulose test
		//String chemicalName = "L-altro-Octo-4,5-diulose";
		




		// need to deal with the following:
		//String chemicalName = "N-[3-[2-[4-(2-methylquinolin-5-yl)piperazin-1-yl]ethyl]phenyl]pyrazine-2-carboxamide";

		//String chemicalName = "N'-methyl-4,7,10,13-tetraoxahexadecanediamide";

		try {
			System.out.println("Parsing chemical name: " + chemicalName);
			System.out.println("==========================================");
			
			// Get NameToStructure instance
			NameToStructure n2s = NameToStructure.getInstance();
			
			// Parse the chemical name and get SMILES
			OpsinResult result = n2s.parseChemicalName(chemicalName);
			
			// Display results
			System.out.println("\nParse Status: " + result.getStatus());
			System.out.println("Output Status: " + result.getOutputStatus());
			
			// Display warnings if any
			if (!result.getWarnings().isEmpty()) {
				System.out.println("\nParse Warnings:");
				for (OpsinWarning warning : result.getWarnings()) {
					System.out.println("  - " + warning.getType() + ": " + warning.getMessage());
				}
			}
			
			if (!result.getOutputWarnings().isEmpty()) {
				System.out.println("\nOutput Warnings:");
				for (OpsinWarning warning : result.getOutputWarnings()) {
					System.out.println("  - " + warning.getType() + ": " + warning.getMessage());
				}
			}
			
			if (result.getStatus() == OpsinResult.OPSIN_RESULT_STATUS.SUCCESS || 
				result.getStatus() == OpsinResult.OPSIN_RESULT_STATUS.WARNING) {
				String smiles = result.getSmiles();
				System.out.println("\nParse SMILES: " + smiles);
				
				// Get output SMILES if available
				if (result.getOutputStructure() != null) {
					try {
						String outputSmiles = SMILESWriter.generateSmiles(result.getOutputStructure(), SmilesOptions.DEFAULT);
						System.out.println("Output SMILES: " + outputSmiles);
					} catch (Exception e) {
						System.out.println("Output SMILES generation failed: " + e.getMessage());
					}
				} else {
					System.out.println("Output SMILES: Not available (output processing failed or not performed)");
					if (result.getOutputStatus() == OpsinResult.OPSIN_RESULT_STATUS.FAILURE) {
						System.out.println("  Reason: " + result.getOutputReasonForFailure());
					}
				}
				
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
		} catch (Exception e) {
			System.err.println("\nERROR OCCURRED:");
			System.err.println("==========================================");
			System.err.println("Error Type: " + e.getClass().getName());
			System.err.println("Error Message: " + e.getMessage());
			System.err.println("\nStack Trace (showing file and line numbers):");
			System.err.println("----------------------------------------");
			StackTraceElement[] stackTrace = e.getStackTrace();
			for (StackTraceElement element : stackTrace) {
				System.err.println("  at " + element.getClassName() + "." + element.getMethodName() + 
					"(" + element.getFileName() + ":" + element.getLineNumber() + ")");
			}
			System.err.println("----------------------------------------");
			e.printStackTrace();
		}
	}
}

