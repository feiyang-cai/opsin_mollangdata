package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Main class to process all test resource files with OPSIN.
 * 
 * This script:
 * 1. Loads all .txt files from the resources directory
 * 2. Reads each file line by line (extracting chemical names)
 * 3. Calls the core function for each chemical name
 * 4. Creates subfolders for each file
 * 5. Saves output logs for each sample
 * 
 * To compile:
 *   javac -cp "target/classes:target/dependency/*" TestOpsinResources.java
 * 
 * To run:
 *   java -cp ".:target/classes:target/dependency/*" uk.ac.cam.ch.wwmm.opsin.TestOpsinResources
 */
public class TestOpsinResources {
	
	/**
	 * Extract chemical name from a line.
	 * Lines may have format: "chemical_name\tInChI=..." or just "chemical_name"
	 * Lines starting with # are comments and should be skipped.
	 * 
	 * @param line Input line from file
	 * @return Chemical name, or null if line should be skipped
	 */
	private static String extractChemicalName(String line) {
		line = line.trim();
		
		// Skip empty lines and comments
		if (line.isEmpty() || line.startsWith("#")) {
			return null;
		}
		
		// Extract chemical name (before tab if present)
		String chemicalName;
		if (line.contains("\t")) {
			chemicalName = line.split("\t")[0].trim();
		} else {
			chemicalName = line.trim();
		}
		
		if (chemicalName.isEmpty()) {
			return null;
		}
		
		return chemicalName;
	}
	
	/**
	 * Create a safe filename from a chemical name.
	 * Replaces problematic characters with underscores.
	 * 
	 * @param chemicalName The chemical name
	 * @return Safe filename
	 */
	private static String createSafeFilename(String chemicalName) {
		StringBuilder sb = new StringBuilder();
		for (char c : chemicalName.toCharArray()) {
			if (Character.isLetterOrDigit(c) || c == '-' || c == '_' || c == '.') {
				sb.append(c);
			} else {
				sb.append('_');
			}
		}
		String result = sb.toString();
		// Limit length
		if (result.length() > 100) {
			result = result.substring(0, 100);
		}
		return result;
	}
	
	/**
	 * Process a single resource file.
	 * Creates a subfolder for the file and processes each chemical name.
	 * 
	 * @param filePath Path to the resource file
	 * @param outputBaseDir Base directory for output logs
	 * @return Object containing statistics and list of failed sample names
	 */
	private static class FileProcessResult {
		int[] stats; // [total, passed, failed]
		java.util.List<String> failedSamples;
		
		FileProcessResult(int[] stats, java.util.List<String> failedSamples) {
			this.stats = stats;
			this.failedSamples = failedSamples;
		}
	}
	
	private static FileProcessResult processResourceFile(File filePath, File outputBaseDir) {
		// Create subfolder for this file
		String fileStem = filePath.getName().replaceFirst("[.][^.]+$", ""); // filename without extension
		File outputDir = new File(outputBaseDir, fileStem);
		outputDir.mkdirs();
		
		System.out.println("Processing file: " + filePath.getName());
		System.out.println("Output directory: " + outputDir.getAbsolutePath());
		
		int sampleCount = 0;
		int passedCount = 0;
		int failedCount = 0;
		java.util.List<String> failedSamples = new java.util.ArrayList<>();
		
		try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
			String line;
			while ((line = reader.readLine()) != null) {
				String chemicalName = extractChemicalName(line);
				
				if (chemicalName == null) {
					continue;
				}
				
				sampleCount++;
				
				// Process the chemical name
				System.out.print("  Processing sample " + sampleCount + ": ");
				if (chemicalName.length() > 50) {
					System.out.print(chemicalName.substring(0, 50) + "...");
				} else {
					System.out.print(chemicalName);
				}
				System.out.println();
				
				ProcessChemicalName.ProcessResult result = ProcessChemicalName.process(chemicalName);
				
				// Determine pass/fail
				if (result.passed) {
					passedCount++;
				} else {
					failedCount++;
					failedSamples.add(chemicalName);
				}
				
				// Create safe filename for output
				String safeName = createSafeFilename(chemicalName);
				String outputFileName = String.format("%04d_%s.txt", sampleCount, safeName);
				File outputFile = new File(outputDir, outputFileName);
				
				// Write output
				ProcessChemicalName.writeOutputToFile(result, outputFile);
			}
		} catch (IOException e) {
			System.err.println("Error reading file " + filePath.getName() + ": " + e.getMessage());
			e.printStackTrace();
			return new FileProcessResult(new int[]{sampleCount, passedCount, failedCount}, failedSamples);
		}
		
		// Create summary file
		File summaryFile = new File(outputDir, "summary.txt");
		try (PrintWriter writer = new PrintWriter(new FileWriter(summaryFile))) {
			writer.println("File: " + filePath.getName());
			writer.println("Total samples: " + sampleCount);
			writer.println("Passed: " + passedCount);
			writer.println("Failed: " + failedCount);
			if (sampleCount > 0) {
				writer.printf("Pass rate: %.2f%%\n", (passedCount * 100.0 / sampleCount));
			} else {
				writer.println("Pass rate: N/A");
			}
		} catch (IOException e) {
			System.err.println("Error writing summary file: " + e.getMessage());
		}
		
		System.out.println("  Completed: " + sampleCount + " samples (" + passedCount + " passed, " + failedCount + " failed)");
		System.out.println();
		
		return new FileProcessResult(new int[]{sampleCount, passedCount, failedCount}, failedSamples);
	}
	
	/**
	 * Find the resources directory.
	 * 
	 * @return Path to resources directory
	 */
	private static File findResourcesDir() {
		// Try to find resources directory relative to current working directory
		// Look for: opsin-core/src/test/resources/uk/ac/cam/ch/wwmm/opsin
		File currentDir = new File(System.getProperty("user.dir"));
		
		// Try relative to current directory
		File resourcesDir = new File(currentDir, "src/test/resources/uk/ac/cam/ch/wwmm/opsin");
		if (resourcesDir.exists() && resourcesDir.isDirectory()) {
			return resourcesDir;
		}
		
		// Try from opsin-core directory
		resourcesDir = new File(currentDir, "opsin-core/src/test/resources/uk/ac/cam/ch/wwmm/opsin");
		if (resourcesDir.exists() && resourcesDir.isDirectory()) {
			return resourcesDir;
		}
		
		// Try absolute path (if we're in opsin-core)
		if (currentDir.getName().equals("opsin-core")) {
			resourcesDir = new File(currentDir, "src/test/resources/uk/ac/cam/ch/wwmm/opsin");
			if (resourcesDir.exists() && resourcesDir.isDirectory()) {
				return resourcesDir;
			}
		}
		
		// Default: assume we're in opsin-core directory
		return new File("src/test/resources/uk/ac/cam/ch/wwmm/opsin");
	}
	
	/**
	 * Main method.
	 * 
	 * @param args Command line arguments (not used)
	 */
	public static void main(String[] args) {
		// Find resources directory
		File resourcesDir = findResourcesDir();
		
		if (!resourcesDir.exists() || !resourcesDir.isDirectory()) {
			System.err.println("Error: Resources directory not found: " + resourcesDir.getAbsolutePath());
			System.err.println("Current working directory: " + System.getProperty("user.dir"));
			System.exit(1);
		}
		
		System.out.println("Resources directory: " + resourcesDir.getAbsolutePath());
		
		// Create output directory
		File outputBaseDir = new File("test_outputs");
		outputBaseDir.mkdirs();
		System.out.println("Output directory: " + outputBaseDir.getAbsolutePath());
		System.out.println();
		
		// Find all .txt files in resources directory
		File[] txtFiles = resourcesDir.listFiles((dir, name) -> name.endsWith(".txt"));
		
		if (txtFiles == null || txtFiles.length == 0) {
			System.err.println("Error: No .txt files found in " + resourcesDir.getAbsolutePath());
			System.exit(1);
		}
		
		// Sort files by name
		java.util.Arrays.sort(txtFiles);
		
		System.out.println("Found " + txtFiles.length + " resource file(s)");
		System.out.println("============================================================");
		System.out.println();
		
		// Process each file
		int totalSamples = 0;
		int totalPassed = 0;
		int totalFailed = 0;
		java.util.Map<String, java.util.List<String>> filesWithFailures = new java.util.LinkedHashMap<>();
		
		for (File txtFile : txtFiles) {
			try {
				FileProcessResult result = processResourceFile(txtFile, outputBaseDir);
				int[] stats = result.stats;
				totalSamples += stats[0];
				totalPassed += stats[1];
				totalFailed += stats[2];
				
				// Track files with failures and their failed sample names
				if (stats[2] > 0) {
					String fileStem = txtFile.getName().replaceFirst("[.][^.]+$", ""); // filename without extension
					filesWithFailures.put(fileStem, result.failedSamples);
				}
			} catch (Exception e) {
				System.err.println("Error processing " + txtFile.getName() + ": " + e.getMessage());
				e.printStackTrace();
				continue;
			}
		}
		
		// Create overall summary
		File overallSummary = new File(outputBaseDir, "overall_summary.txt");
		try (PrintWriter writer = new PrintWriter(new FileWriter(overallSummary))) {
			writer.println("Overall Test Summary");
			writer.println("============================================================");
			writer.println("Total files processed: " + txtFiles.length);
			writer.println("Total samples: " + totalSamples);
			writer.println("Total passed: " + totalPassed);
			writer.println("Total failed: " + totalFailed);
			if (totalSamples > 0) {
				writer.printf("Overall pass rate: %.2f%%\n", (totalPassed * 100.0 / totalSamples));
			} else {
				writer.println("Overall pass rate: N/A");
			}
			writer.println();
			writer.println("File directories with failures:");
			if (filesWithFailures.isEmpty()) {
				writer.println("  (None - all tests passed)");
			} else {
				for (java.util.Map.Entry<String, java.util.List<String>> entry : filesWithFailures.entrySet()) {
					String fileDir = entry.getKey();
					java.util.List<String> failedSamples = entry.getValue();
					writer.println("  - " + fileDir);
					for (String sampleName : failedSamples) {
						writer.println("    - " + sampleName);
					}
				}
			}
		} catch (IOException e) {
			System.err.println("Error writing overall summary: " + e.getMessage());
		}
		
		System.out.println("============================================================");
		System.out.println("Processing complete!");
		System.out.println("Total samples: " + totalSamples);
		System.out.println("Passed: " + totalPassed);
		System.out.println("Failed: " + totalFailed);
		if (totalSamples > 0) {
			System.out.printf("Pass rate: %.2f%%\n", (totalPassed * 100.0 / totalSamples));
		}
		System.out.println();
		System.out.println("Results saved to: " + outputBaseDir.getAbsolutePath());
		System.out.println("Overall summary: " + overallSummary.getAbsolutePath());
	}
}

