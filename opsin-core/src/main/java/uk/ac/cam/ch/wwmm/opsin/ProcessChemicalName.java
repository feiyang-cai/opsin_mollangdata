package uk.ac.cam.ch.wwmm.opsin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Core class to process a single chemical name with OPSIN.
 * This class provides a reusable method that can be called from other classes.
 * 
 * Similar to ExampleGetSmilesAndXml but returns structured results.
 */
public class ProcessChemicalName {
	
	/**
	 * Result class to hold all information about processing a chemical name.
	 */
	public static class ProcessResult {
		public String chemicalName;
		public boolean passed;
		public String reason;
		public String parseStatus;
		public String outputStatus;
		public String parseSmiles;
		public String outputSmiles;
		public String parseXml;
		public String outputXml;
		public List<String> parseWarnings = new ArrayList<>();
		public List<String> outputWarnings = new ArrayList<>();
		public String errorMessage;
		public String parseMessage;  // Message from getMessage() - includes warnings or failure reason
		public String outputMessage;  // Message from getOutputMessage() - includes output warnings
		public String outputReasonForFailure;  // Specific reason for output failure
	}
	
	/**
	 * Process a single chemical name using OPSIN.
	 * 
	 * This method mimics the behavior of ExampleGetSmilesAndXml:
	 * - Gets parse SMILES and output SMILES
	 * - Gets parse XML and output XML
	 * - Compares SMILES
	 * - Checks status (SUCCESS, WARNING, or FAILURE)
	 * - Returns all information including pass/fail status
	 * 
	 * @param chemicalName The IUPAC chemical name to process
	 * @return ProcessResult containing all information
	 */
	public static ProcessResult process(String chemicalName) {
		ProcessResult result = new ProcessResult();
		result.chemicalName = chemicalName;
		result.passed = false;
		
		try {
			// Get NameToStructure instance (like ExampleGetSmilesAndXml)
			NameToStructure n2s = NameToStructure.getInstance();
			
			// Parse the chemical name
			OpsinResult opsinResult = n2s.parseChemicalName(chemicalName);
			
			// Extract status
			result.parseStatus = opsinResult.getStatus().toString();
			result.outputStatus = opsinResult.getOutputStatus().toString();
			
			// Extract messages (these include warnings or failure reasons)
			result.parseMessage = opsinResult.getMessage();
			result.outputMessage = opsinResult.getOutputMessage();
			result.outputReasonForFailure = opsinResult.getOutputReasonForFailure();
			
			// Extract warnings
			for (OpsinWarning warning : opsinResult.getWarnings()) {
				result.parseWarnings.add(warning.getType() + ": " + warning.getMessage());
			}
			
			for (OpsinWarning warning : opsinResult.getOutputWarnings()) {
				result.outputWarnings.add(warning.getType() + ": " + warning.getMessage());
			}
			
			// Check if parsing was successful or had warnings
			if (opsinResult.getStatus() == OpsinResult.OPSIN_RESULT_STATUS.SUCCESS || 
				opsinResult.getStatus() == OpsinResult.OPSIN_RESULT_STATUS.WARNING) {
				
				// Get parse SMILES
				result.parseSmiles = opsinResult.getSmiles();
				
				// Get parse XML
				result.parseXml = opsinResult.getParseXml();
				
				// Get output SMILES if available
				if (opsinResult.getOutputStructure() != null) {
					try {
						result.outputSmiles = SMILESWriter.generateSmiles(
							opsinResult.getOutputStructure(), 
							SmilesOptions.DEFAULT
						);
					} catch (Exception e) {
						result.outputSmiles = null;
						result.errorMessage = "Output SMILES generation failed: " + e.getMessage();
					}
				}
				
				// Get output XML
				result.outputXml = opsinResult.getOutputXml();
				
				// Determine if passed
				// Pass if: both statuses are SUCCESS AND SMILES match
				if (result.parseStatus.equals("SUCCESS") && result.outputStatus.equals("SUCCESS")) {
					if (result.parseSmiles != null && result.outputSmiles != null) {
						if (result.parseSmiles.equals(result.outputSmiles)) {
							result.passed = true;
							result.reason = "PASS";
						} else {
							result.passed = false;
							result.reason = "SMILES mismatch: parse='" + result.parseSmiles + 
								"' vs output='" + result.outputSmiles + "'";
						}
					} else if (result.parseSmiles != null && result.outputSmiles == null) {
						result.passed = false;
						result.reason = "Output SMILES not available";
					} else {
						result.passed = false;
						result.reason = "Parse SMILES not available";
					}
				} else {
					// Status is not SUCCESS (could be WARNING or FAILURE)
					StringBuilder statusReason = new StringBuilder();
					if (!result.parseStatus.equals("SUCCESS")) {
						statusReason.append("Parse status: ").append(result.parseStatus);
						// Include the actual message/warning/error
						if (result.parseMessage != null && !result.parseMessage.isEmpty()) {
							statusReason.append(" (").append(result.parseMessage).append(")");
						}
					}
					if (!result.outputStatus.equals("SUCCESS")) {
						if (statusReason.length() > 0) {
							statusReason.append("; ");
						}
						statusReason.append("Output status: ").append(result.outputStatus);
						// Include the actual message/warning/error
						if (result.outputStatus.equals("FAILURE")) {
							// For FAILURE, use the specific failure reason
							if (result.outputReasonForFailure != null && !result.outputReasonForFailure.isEmpty()) {
								statusReason.append(" (").append(result.outputReasonForFailure).append(")");
							} else if (result.outputMessage != null && !result.outputMessage.isEmpty()) {
								statusReason.append(" (").append(result.outputMessage).append(")");
							}
						} else if (result.outputMessage != null && !result.outputMessage.isEmpty()) {
							// For WARNING, use the output message
							statusReason.append(" (").append(result.outputMessage).append(")");
						}
					}
					result.passed = false;
					result.reason = statusReason.toString();
				}
			} else {
				// Parsing failed
				result.errorMessage = opsinResult.getMessage();
				result.parseMessage = opsinResult.getMessage();
				result.passed = false;
				result.reason = "Parse failed: " + opsinResult.getMessage();
			}
		} catch (Exception e) {
			result.parseStatus = "ERROR";
			result.outputStatus = "ERROR";
			result.errorMessage = e.getMessage();
			result.passed = false;
			result.reason = "Exception occurred: " + e.getMessage();
		}
		
		return result;
	}
	
	/**
	 * Format the result into a readable output string.
	 * Output format matches ExampleGetSmilesAndXml but with output XML first.
	 * 
	 * @param result The ProcessResult to format
	 * @return Formatted string with all information
	 */
	public static String formatOutput(ProcessResult result) {
		StringBuilder sb = new StringBuilder();
		
		// First line: PASS or NOT PASS with reason
		if (result.passed) {
			sb.append("PASS\n");
		} else {
			sb.append("NOT PASS: ").append(result.reason).append("\n");
		}
		
		sb.append("\n");
		sb.append("Chemical Name: ").append(result.chemicalName).append("\n");
		sb.append("==================================================\n");
		sb.append("\n");
		sb.append("Parse Status: ").append(result.parseStatus).append("\n");
		sb.append("Output Status: ").append(result.outputStatus).append("\n");
		
		// Warnings
		if (!result.parseWarnings.isEmpty()) {
			sb.append("\n");
			sb.append("Parse Warnings:\n");
			for (String warning : result.parseWarnings) {
				sb.append("  - ").append(warning).append("\n");
			}
		}
		
		if (!result.outputWarnings.isEmpty()) {
			sb.append("\n");
			sb.append("Output Warnings:\n");
			for (String warning : result.outputWarnings) {
				sb.append("  - ").append(warning).append("\n");
			}
		}
		
		// SMILES
		sb.append("\n");
		if (result.parseSmiles != null) {
			sb.append("Parse SMILES: ").append(result.parseSmiles).append("\n");
		} else {
			sb.append("Parse SMILES: Not available\n");
		}
		
		if (result.outputSmiles != null) {
			sb.append("Output SMILES: ").append(result.outputSmiles).append("\n");
		} else {
			sb.append("Output SMILES: Not available\n");
			if (result.outputStatus.equals("FAILURE")) {
				// Show the specific failure reason if available
				String failureReason = result.outputReasonForFailure != null && !result.outputReasonForFailure.isEmpty() 
					? result.outputReasonForFailure 
					: (result.outputMessage != null && !result.outputMessage.isEmpty() 
						? result.outputMessage 
						: (result.errorMessage != null ? result.errorMessage : "Unknown"));
				sb.append("  Reason: ").append(failureReason).append("\n");
			}
		}
		
		// Show parse message if status is WARNING or FAILURE
		if (!result.parseStatus.equals("SUCCESS") && result.parseMessage != null && !result.parseMessage.isEmpty()) {
			sb.append("\n");
			sb.append("Parse Message: ").append(result.parseMessage).append("\n");
		}
		
		// Show output message if status is WARNING or FAILURE
		if (!result.outputStatus.equals("SUCCESS") && result.outputMessage != null && !result.outputMessage.isEmpty()) {
			sb.append("\n");
			sb.append("Output Message: ").append(result.outputMessage).append("\n");
		}
		
		// Show output failure reason if status is FAILURE
		if (result.outputStatus.equals("FAILURE") && result.outputReasonForFailure != null && !result.outputReasonForFailure.isEmpty()) {
			sb.append("\n");
			sb.append("Output Failure Reason: ").append(result.outputReasonForFailure).append("\n");
		}
		
		// XML (output XML first as requested)
		if (result.outputXml != null) {
			sb.append("\n");
			sb.append("Output XML:\n");
			sb.append("----------------------------------------\n");
			sb.append(result.outputXml).append("\n");
			sb.append("----------------------------------------\n");
		}
		
		if (result.parseXml != null) {
			sb.append("\n");
			sb.append("Parse XML:\n");
			sb.append("----------------------------------------\n");
			sb.append(result.parseXml).append("\n");
			sb.append("----------------------------------------\n");
		}
		
		// Error message if any
		if (result.errorMessage != null && !result.errorMessage.isEmpty()) {
			sb.append("\n");
			sb.append("Error: ").append(result.errorMessage).append("\n");
		}
		
		return sb.toString();
	}
	
	/**
	 * Write the formatted output to a file.
	 * 
	 * @param result The ProcessResult to write
	 * @param outputFile The file to write to
	 * @throws IOException If writing fails
	 */
	public static void writeOutputToFile(ProcessResult result, File outputFile) throws IOException {
		try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile))) {
			writer.print(formatOutput(result));
		}
	}
}

