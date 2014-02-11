package owl.core.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * Some utilities to query jobs from the Sun Grid Engine system.
 * The java API DRMAA will do most of it but there are some 
 * limitations in DRMAA, for instance querying for jobs that are
 * already finished, which can only be done with qacct.
 * 
 * @author duarte_j
 *
 */
public class GEUtil {
	
	private static final Pattern FAILED_REGEX = Pattern.compile("^failed\\s+(\\d+)\\s+.*$");
	private static final Pattern EXITSTATUS_REGEX = Pattern.compile("^exit_status\\s+(\\d+)\\s+.*$");
	
	/**
	 * Given a Grid Engine job id or job name returns the GEExitStatus by using qacct
	 * GEExitStatus is composed of 2 parts: 
	 * - failed: GE's internal indicator that a job runs or fails to run (for instance before submission)
	 * - exit_status: the actual exit status of the job itself
	 * @param qacctExec
	 * @param job the GE job id or the job name (passed with -N, must start with a letter)
	 * @return the GEExitStatus or null if given job id/name does not exist
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static GEExitStatus getExitStatus(File qacctExec, String job) throws IOException, InterruptedException {
		
		String[] cmd = {qacctExec.getAbsolutePath(), "-j", job};
		
		Process proc = Runtime.getRuntime().exec(cmd);

		int failed = -1;
		int exitStatus = -1;
		
		BufferedReader outBR = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		String line;
		while((line = outBR.readLine()) != null) {
			if (line.startsWith("failed")) {
				Matcher m = FAILED_REGEX.matcher(line);
				if (m.matches()) {
					failed = Integer.parseInt(m.group(1));
				}
				continue;
			}
			if (line.startsWith("exit_status")) {
				Matcher m = EXITSTATUS_REGEX.matcher(line);
				if (m.matches()) {
					exitStatus = Integer.parseInt(m.group(1));
				}				
				break;
			}
		}
		boolean unknownJob = false;
		BufferedReader errBR = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
		while((line = errBR.readLine()) != null) {
			if (line.matches("^error: job .* not found$")) {
				unknownJob = true;
				break;
			}
		}

		int exitValue = proc.waitFor();
 
		if (exitValue!=0) {
			if (unknownJob) {
				return null;
			} else {
				throw new IOException("Exit value of qacct was "+exitValue);
			}
		}
		
		return new GEExitStatus(failed, exitStatus); 
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		for (int i=1;i<100;i++) {
			GEExitStatus es = getExitStatus(new File("/usr/bin/qacct"),String.format("%d", i));
		
			if (es==null) {
				System.out.println("No such job id/name "+i);
				continue;
			}
			System.out.println("Failed: "+es.getFailed());
			System.out.println("Exit status: "+es.getExitStatus());
		}
	}
}
