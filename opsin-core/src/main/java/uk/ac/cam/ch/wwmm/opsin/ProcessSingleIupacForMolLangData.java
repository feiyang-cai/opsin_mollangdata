package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Utility to process ONE IUPAC name and output:
 * - OPSIN SMILES/XML ("parse")
 * - MolLangData SMILES/XML ("output structure")
 * - statuses + warnings/failure reasons for both
 * - final status
 *
 * Final PASS requires:
 *  1) OPSIN status == SUCCESS
 *  2) MolLangData status == SUCCESS
 *  3) canonical(opsinSmiles) == canonical(mollangdataSmiles)
 *
 * Notes:
 * - SMILES canonicalisation is done using OPSIN's own SMILES parser + SMILESWriter, so this
 *   comparison is robust to different SMILES atom orderings.
 *
 * IMPORTANT:
 * - This class does NOT write any files. It prints a single JSON object to stdout.
 *   Your Python/CLI wrapper can decide what to persist.
 */
public class ProcessSingleIupacForMolLangData {

	private static class EvalResult {
		ProcessChemicalName.ProcessResult opsin;

		String canonicalParseSmiles;
		String canonicalOutputSmiles;

		boolean parseAndOutputSmilesMatch;

		boolean finalPass;
		String finalStatus;
	}

	private static String canonicalizeSmilesOrNull(String smiles, boolean[] parsedByOpsinOut) {
		if (parsedByOpsinOut != null && parsedByOpsinOut.length > 0) {
			parsedByOpsinOut[0] = false;
		}
		if (smiles == null) {
			return null;
		}
		String trimmed = smiles.trim();
		if (trimmed.isEmpty()) {
			return null;
		}
		try {
			Fragment frag = new SMILESFragmentBuilder(new IDManager()).build(trimmed);
			if (parsedByOpsinOut != null && parsedByOpsinOut.length > 0) {
				parsedByOpsinOut[0] = true;
			}
			return SMILESWriter.generateSmiles(frag, SmilesOptions.DEFAULT);
		}
		catch (Exception e) {
			// Fall back to raw string; still useful for debugging, but comparison may be less reliable.
			return trimmed;
		}
	}

	private static boolean equalNullable(String a, String b) {
		if (a == null && b == null) {
			return true;
		}
		if (a == null || b == null) {
			return false;
		}
		return a.equals(b);
	}

	private static EvalResult evaluate(String iupac) {
		EvalResult er = new EvalResult();
		er.opsin = ProcessChemicalName.process(iupac);

		er.canonicalParseSmiles = canonicalizeSmilesOrNull(er.opsin.parseSmiles, null);
		er.canonicalOutputSmiles = canonicalizeSmilesOrNull(er.opsin.outputSmiles, null);

		er.parseAndOutputSmilesMatch = equalNullable(er.canonicalParseSmiles, er.canonicalOutputSmiles);

		boolean parseSuccess = "SUCCESS".equals(er.opsin.parseStatus);
		boolean outputSuccess = "SUCCESS".equals(er.opsin.outputStatus);

		// Final status rules:
		// 1) PASS: both SUCCESS and SMILES match
		// 2) Opsin warning/fail: parse is WARNING/FAIL (anything not SUCCESS is acceptable here)
		// 3) MolLangData warning/fail: parse SUCCESS but output WARNING/FAIL
		er.finalPass = parseSuccess && outputSuccess && er.parseAndOutputSmilesMatch;
		if (er.finalPass) {
			er.finalStatus = "PASS";
			return er;
		}
		if (!parseSuccess) {
			er.finalStatus = "OPSIN warning/fail";
			return er;
		}
		if (!outputSuccess) {
			er.finalStatus = "MolLangData warning/fail";
			return er;
		}
		// Both SUCCESS but SMILES mismatch
		er.finalStatus = "MolLangData SMILES mismatch";
		return er;
	}

	static Map<String, Object> processOneToJson(String iupac) {
		EvalResult er = evaluate(iupac);

		Map<String, Object> json = new LinkedHashMap<>();
		json.put("iupac", iupac);

		// Status + messages
		json.put("opsin_status", er.opsin.parseStatus);
		json.put("mollangdata_status", er.opsin.outputStatus);
		json.put("opsin_message", er.opsin.parseMessage);
		json.put("mollangdata_message", er.opsin.outputMessage);
		json.put("mollangdata_reason_for_failure", er.opsin.outputReasonForFailure);

		// Warnings
		json.put("opsin_warnings", er.opsin.parseWarnings == null ? new ArrayList<>() : er.opsin.parseWarnings);
		json.put("mollangdata_warnings", er.opsin.outputWarnings == null ? new ArrayList<>() : er.opsin.outputWarnings);

		// Artifacts (as strings; caller can write to files if desired)
		json.put("opsin_smiles", er.opsin.parseSmiles);
		json.put("mollangdata_smiles", er.opsin.outputSmiles);
		json.put("opsin_xml", er.opsin.parseXml);
		json.put("mollangdata_xml", er.opsin.outputXml);

		// Canonical comparisons
		json.put("opsin_smiles_canonical", er.canonicalParseSmiles);
		json.put("mollangdata_smiles_canonical", er.canonicalOutputSmiles);
		json.put("opsin_mollangdata_smiles_match", er.parseAndOutputSmilesMatch);

		// Final status
		json.put("final_pass", er.finalPass);
		json.put("final_status", er.finalStatus);

		return json;
	}

	/**
	 * Minimal JSON writer (no external deps).
	 */
	static String toJson(Map<String, Object> m) {
		StringBuilder sb = new StringBuilder();
		sb.append("{");
		boolean first = true;
		for (Map.Entry<String, Object> e : m.entrySet()) {
			if (!first) sb.append(",");
			first = false;
			sb.append("\"").append(escapeJson(e.getKey())).append("\":");
			sb.append(valueToJson(e.getValue()));
		}
		sb.append("}");
		return sb.toString();
	}

	private static String valueToJson(Object v) {
		if (v == null) return "null";
		if (v instanceof Boolean) return ((Boolean) v) ? "true" : "false";
		if (v instanceof Number) return v.toString();
		if (v instanceof Map) {
			@SuppressWarnings("unchecked")
			Map<String, Object> m = (Map<String, Object>) v;
			return toJson(m);
		}
		if (v instanceof Iterable) {
			StringBuilder sb = new StringBuilder();
			sb.append("[");
			boolean first = true;
			for (Object x : (Iterable<?>) v) {
				if (!first) sb.append(",");
				first = false;
				sb.append(valueToJson(x));
			}
			sb.append("]");
			return sb.toString();
		}
		return "\"" + escapeJson(v.toString()) + "\"";
	}

	private static String escapeJson(String s) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			switch (c) {
				case '\\': sb.append("\\\\"); break;
				case '"': sb.append("\\\""); break;
				case '\n': sb.append("\\n"); break;
				case '\r': sb.append("\\r"); break;
				case '\t': sb.append("\\t"); break;
				default:
					if (c < 0x20) sb.append(String.format("\\u%04x", (int) c));
					else sb.append(c);
			}
		}
		return sb.toString();
	}

	private static void usage() {
		System.err.println("Usage:");
		System.err.println("  java ... uk.ac.cam.ch.wwmm.opsin.ProcessSingleIupacForMolLangData \"<IUPAC>\"");
		System.err.println();
		System.err.println("Output:");
		System.err.println("  - Prints a single JSON object to stdout (no files written).");
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 1) {
			usage();
			System.exit(1);
		}

		try {
			String iupac = args[0];

			Map<String, Object> json = processOneToJson(iupac);

			System.out.print(toJson(json));
			System.out.print("\n");
		}
		catch (Exception e) {
			Map<String, Object> err = new LinkedHashMap<>();
			err.put("error", true);
			err.put("exception_class", e.getClass().getName());
			err.put("exception_message", e.getMessage());
			System.out.print(toJson(err));
			System.out.print("\n");
			System.exit(2);
		}
	}
}


