package gps.test;

import gps.acquisition.Acquisition;

public class AcquisitionTest {

	public static void main(String[] args) {
		

		Acquisition acq = null;
		AcquisitionTestCase testCase = new AcquisitionTestCase1();
		
		boolean acquisition = false;
		int codeVersch = -1;
		

		int nrOfSamples = 2000;
		acq = new Acquisition(nrOfSamples);
		
		float[][] inputSamples = testCase.getInputSamples();
		for(int i = 0; i < nrOfSamples; i++){
			float real = inputSamples[i][0];
			float imag = inputSamples[i][1];
			acq.enterSample(real, imag);
		}
		
		float[] inputCodes = testCase.getInputCodes();
		for(int i = 0; i < nrOfSamples; i++){
			float real = inputCodes[2*i];
			float imag = inputCodes[2*i+1];
			acq.enterCode(real, imag);
		}

		acquisition = testCase.getAcquisitionMessage();
		codeVersch = testCase.getCodePhase();

		boolean res = acq.startAcquisition(testCase.getDopplerFrequency());
		
		boolean passed = res == acquisition && acq.getCodePhase()== codeVersch;
		
		System.out.println((passed?"PASSED":"FAILED"));
		if(!passed){
			System.out.println("Epected " + acquisition + " acquistion");
			System.out.println("    " + codeVersch);
			
			System.out.println("Got " + res + " acquistion");
			System.out.println("     " + acq.getCodePhase());
		}

	}
	

}
