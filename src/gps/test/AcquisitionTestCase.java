package gps.test;

public interface AcquisitionTestCase {
	
	public float[] getInputCodes();
	
	public float[][] getInputSamples();
	
	public int getDopplerFrequency();
	public boolean getAcquisitionMessage();
	public int getCodePhase();
	public float getMaxMagnitude();
	public float getInputPower();
	public float getGamma();
}
