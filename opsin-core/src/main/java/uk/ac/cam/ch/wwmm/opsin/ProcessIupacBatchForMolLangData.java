package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Batch processor for many IUPAC names.
 *
 * Input: UTF-8 lines on stdin, one IUPAC per line.
 * Output: one JSON object per input line (JSONL), written to stdout.
 *
 * This avoids JVM startup overhead per molecule and is intended for Python wrappers.
 */
public class ProcessIupacBatchForMolLangData {

	public static void main(String[] args) throws Exception {
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in, StandardCharsets.UTF_8));
		PrintWriter out = new PrintWriter(System.out, true);

		String line;
		while ((line = br.readLine()) != null) {
			String iupac = line.trim();
			if (iupac.isEmpty()) {
				continue;
			}
			try {
				Map<String, Object> json = ProcessSingleIupacForMolLangData.processOneToJson(iupac);
				out.println(ProcessSingleIupacForMolLangData.toJson(json));
			}
			catch (Exception e) {
				Map<String, Object> err = new LinkedHashMap<>();
				err.put("error", true);
				err.put("iupac", iupac);
				err.put("exception_class", e.getClass().getName());
				err.put("exception_message", e.getMessage());
				out.println(ProcessSingleIupacForMolLangData.toJson(err));
			}
		}
		out.flush();
	}
}

