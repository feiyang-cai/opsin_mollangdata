package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;

/**
 * Generates SMILES and XML outputs for all test cases in the test resources.
 * For each test case file, creates a subfolder and saves each test case's output
 * in a separate file within that subfolder.
 */
public class GenerateTestOutputs {
	
	private static final String RESOURCES_DIR = "uk/ac/cam/ch/wwmm/opsin/";
	private static final String OUTPUT_BASE_DIR = "test_outputs/";
	
	public static void main(String[] args) {
		System.out.println("Starting test case output generation...");
		System.out.println("==========================================");
		
		// Initialize OPSIN
		NameToStructure n2s = NameToStructure.getInstance();
		
		// Get the resources directory
		File resourcesDir = null;
		
		// Try to get from classloader first (works when running from compiled classes)
		ClassLoader classLoader = GenerateTestOutputs.class.getClassLoader();
		URL resourcesUrl = classLoader.getResource(RESOURCES_DIR);
		
		if (resourcesUrl != null) {
			try {
				// Handle both file system and JAR resources
				if (!resourcesUrl.getProtocol().equals("jar")) {
					resourcesDir = new File(resourcesUrl.toURI());
				}
			} catch (Exception e) {
				// Will try fallback below
			}
		}
		
		// Fallback: try to read from source directory structure
		if (resourcesDir == null || !resourcesDir.exists()) {
			// Try relative path from project root
			File testResourcesDir = new File("opsin-inchi/src/test/resources/" + RESOURCES_DIR);
			if (testResourcesDir.exists() && testResourcesDir.isDirectory()) {
				resourcesDir = testResourcesDir;
			} else {
				// Try absolute path based on current working directory
				File currentDir = new File(".").getAbsoluteFile();
				File parentDir = currentDir.getParentFile();
				if (parentDir != null) {
					testResourcesDir = new File(parentDir, "opsin-inchi/src/test/resources/" + RESOURCES_DIR);
					if (testResourcesDir.exists() && testResourcesDir.isDirectory()) {
						resourcesDir = testResourcesDir;
					}
				}
			}
		}
		
		if (resourcesDir == null || !resourcesDir.exists()) {
			System.err.println("Could not find resources directory: " + RESOURCES_DIR);
			System.err.println("Tried classloader resource and file system paths.");
			return;
		}
		
		// Create output base directory
		File outputBaseDir = new File(OUTPUT_BASE_DIR);
		if (!outputBaseDir.exists()) {
			outputBaseDir.mkdirs();
		}
		
		// Process all .txt files in the resources directory
		File[] testFiles = resourcesDir.listFiles((dir, name) -> name.endsWith(".txt"));
		
		if (testFiles == null || testFiles.length == 0) {
			System.err.println("No .txt files found in resources directory: " + resourcesDir.getAbsolutePath());
			return;
		}
		
		int totalFiles = testFiles.length;
		int processedFiles = 0;
		int totalTestCases = 0;
		int successfulCases = 0;
		int failedCases = 0;
		
		for (File testFile : testFiles) {
			String fileName = testFile.getName();
			String baseName = fileName.substring(0, fileName.lastIndexOf('.'));
			
			System.out.println("\nProcessing file: " + fileName);
			System.out.println("----------------------------------------");
			
			// Create subfolder for this test file
			File subfolder = new File(outputBaseDir, baseName);
			if (!subfolder.exists()) {
				subfolder.mkdirs();
			}
			
			int caseNumber = 0;
			int fileSuccessCount = 0;
			int fileFailCount = 0;
			
			try (BufferedReader reader = new BufferedReader(new FileReader(testFile))) {
				String line;
				while ((line = reader.readLine()) != null) {
					line = line.trim();
					
					// Skip empty lines and comments
					if (line.isEmpty() || line.startsWith("#")) {
						continue;
					}
					
					// Split by tab to get chemical name (first column)
					String[] parts = line.split("\t");
					if (parts.length == 0 || parts[0].trim().isEmpty()) {
						continue;
					}
					
					String chemicalName = parts[0].trim();
					caseNumber++;
					totalTestCases++;
					
					// Parse the chemical name (using default configuration)
					OpsinResult result = n2s.parseChemicalName(chemicalName);
					
					// Create output file for this test case
					String outputFileName = String.format("%04d_%s.txt", caseNumber, sanitizeFileName(chemicalName));
					File outputFile = new File(subfolder, outputFileName);
					
					try (FileWriter writer = new FileWriter(outputFile)) {
						writer.write("Chemical Name: " + chemicalName + "\n");
						writer.write("==========================================\n\n");
						
						writer.write("Status: " + result.getStatus() + "\n");
						
						if (result.getStatus() == OpsinResult.OPSIN_RESULT_STATUS.SUCCESS) {
							fileSuccessCount++;
							successfulCases++;
							
							// Write SMILES
							String smiles = result.getSmiles();
							if (smiles != null) {
								writer.write("\nSMILES:\n");
								writer.write("----------------------------------------\n");
								writer.write(smiles + "\n");
								writer.write("----------------------------------------\n");
							}
							
							// Write Parse XML
							String parseXml = result.getParseXml();
							if (parseXml != null) {
								writer.write("\nParse XML:\n");
								writer.write("----------------------------------------\n");
								writer.write(parseXml + "\n");
								writer.write("----------------------------------------\n");
							}
							
							// Write Output XML
							String outputXml = result.getOutputXml();
							if (outputXml != null) {
								writer.write("\nOutput XML:\n");
								writer.write("----------------------------------------\n");
								writer.write(outputXml + "\n");
								writer.write("----------------------------------------\n");
							}
							
							// Also write CML if available
							String cml = result.getCml();
							if (cml != null) {
								writer.write("\nCML:\n");
								writer.write("----------------------------------------\n");
								writer.write(cml + "\n");
								writer.write("----------------------------------------\n");
							}
							
						} else if (result.getStatus() == OpsinResult.OPSIN_RESULT_STATUS.WARNING) {
							fileSuccessCount++;
							successfulCases++;
							
							writer.write("\nWarning: " + result.getMessage() + "\n");
							
							String smiles = result.getSmiles();
							if (smiles != null) {
								writer.write("\nSMILES:\n");
								writer.write("----------------------------------------\n");
								writer.write(smiles + "\n");
								writer.write("----------------------------------------\n");
							}
							
							String parseXml = result.getParseXml();
							if (parseXml != null) {
								writer.write("\nParse XML:\n");
								writer.write("----------------------------------------\n");
								writer.write(parseXml + "\n");
								writer.write("----------------------------------------\n");
							}
							
							String outputXml = result.getOutputXml();
							if (outputXml != null) {
								writer.write("\nOutput XML:\n");
								writer.write("----------------------------------------\n");
								writer.write(outputXml + "\n");
								writer.write("----------------------------------------\n");
							}
							
						} else {
							fileFailCount++;
							failedCases++;
							
							writer.write("\nParsing failed: " + result.getMessage() + "\n");
						}
					}
					
					if (caseNumber % 10 == 0) {
						System.out.print(".");
					}
				}
			} catch (IOException e) {
				System.err.println("\nError reading file " + fileName + ": " + e.getMessage());
				continue;
			}
			
			processedFiles++;
			System.out.println("\nCompleted: " + caseNumber + " test cases");
			System.out.println("  Success: " + fileSuccessCount + ", Failed: " + fileFailCount);
		}
		
		System.out.println("\n\n==========================================");
		System.out.println("Summary:");
		System.out.println("  Files processed: " + processedFiles + " / " + totalFiles);
		System.out.println("  Total test cases: " + totalTestCases);
		System.out.println("  Successful: " + successfulCases);
		System.out.println("  Failed: " + failedCases);
		System.out.println("  Output directory: " + new File(OUTPUT_BASE_DIR).getAbsolutePath());
		System.out.println("==========================================");
	}
	
	/**
	 * Sanitizes a string to be used as a filename by removing/replacing invalid characters.
	 */
	private static String sanitizeFileName(String name) {
		// Replace invalid filename characters with underscores
		String sanitized = name.replaceAll("[^a-zA-Z0-9._-]", "_");
		// Limit length to avoid filesystem issues
		if (sanitized.length() > 100) {
			sanitized = sanitized.substring(0, 100);
		}
		return sanitized;
	}
}

