package gps.acquisition;

import cgra.pe.PETrigonometry;

public class Acquisition {

    private int CodeVerschiebung=0;
    // threshold
    private float Gamma_r=0.015f;

    // set within the constructor function
    private int nrOfSamples;
    private float[][] input_data;
    private float[][] input_code;
    private int count=0;
    private float S_max=0.0f;
    private float P_in=0.0f;

    /************************** constructor *************************************/
    public Acquisition(int nrOfSamples) {
        // TODO Auto-generated constructor stub
        this.nrOfSamples = nrOfSamples;
        this.input_data = new float[this.nrOfSamples][2];
        this.input_code = new float[this.nrOfSamples][2];
    }

    /************************** enterSample *************************************/
    public void enterSample(float real, float imag){
        // X_in_real
        input_data[count][0] = real;
        // X_in_imag
        input_data[count][1] = imag;
        count ++;
        if (count == nrOfSamples){
            count = 0;
        }
    }

    /************************** enterCode *************************************/
    public void enterCode(float real, float imag){
        input_code[count][0] = real;
        input_code[count][1] = imag;
        count = count+1;
    }

    /************************** startAcquisition *************************************/
    public boolean startAcquisition(int dopplerFrequency){

        final float PI = 3.14159f;
        float Gamma_R = 0.0f;
        int F_s = 2000000;

        /************************** step1: X_fd calculation *************************************/
        float num_X_fd = 2*PI*dopplerFrequency/F_s;
        float[][] X_fd = new float[nrOfSamples][2];
        X_fd = X_fd_Cal(num_X_fd);
        
        /************************** step2: DFT(X_fd) & DFT(C) calculation *************************************/
        //float multiply[][] = new float[nrOfSamples][2];
        float num_DFT = 2*PI/nrOfSamples;
        float multiply[][] = DFT_Cal(num_DFT, X_fd);
        
        /************************** step4: R_fd calculation *************************************/
        float num_R_fd = 2*PI/nrOfSamples;
        R_fd_Cal(num_R_fd, multiply);
      
        P_in = P_in*nrOfSamples;

        Gamma_R = S_max/P_in;
        if(Gamma_R>Gamma_r){
            return true;
        }else {
            return false;
        }

    }

	private float[][] X_fd_Cal(float num_X_fd) {
		// TODO Auto-generated method stub
		float real_part = 0.0f;
		float imag_part = 0.0f;
		float factor = 0.0f;
		float cosValue =0.0f;
		float sinValue =0.0f;
		float[][] X_fd = new float[nrOfSamples][2];
		for(int i=0;i<nrOfSamples;i=i+1){
            real_part = input_data[i][0];
            imag_part = input_data[i][1];
            P_in = P_in+ real_part*real_part + imag_part*imag_part;
            factor = num_X_fd*i;
            cosValue = PETrigonometry.cos(factor);
            sinValue = PETrigonometry.sin(factor);
            X_fd[i][0] = real_part*cosValue + imag_part*sinValue;
            X_fd[i][1] = imag_part*cosValue - real_part*sinValue;
        }
		return X_fd;
	}

	private float[][] DFT_Cal(float num_DFT, float[][] X_fd) {
		// TODO Auto-generated method stub
		int numSample = nrOfSamples;
		float factor = 0.0f;	
		float[][] DFT_data = new float[numSample][2];
		float[][] DFT_code = new float[numSample][2];
		float multiply[][] = new float[numSample][2];
		for (int i=0; i<numSample; i=i+4) {
            float cosValue =0.0f;
            float sinValue =0.0f;
            for (int j=0; j<numSample; j=j+1) {

                // DFT_X_fd
            	factor = num_DFT*i*j;
            	cosValue = PETrigonometry.cos(factor);
				sinValue = PETrigonometry.sin(factor);	
				DFT_data[i][0] +=  X_fd[j][0] * cosValue + X_fd[j][1] * sinValue;
				DFT_data[i][1] += -X_fd[j][0] * sinValue + X_fd[j][1] * cosValue;
                // DFT_C
				DFT_code[i][0] +=  input_code[j][0] * cosValue + input_code[j][1] * sinValue;
                DFT_code[i][1] += -input_code[j][0] * sinValue + input_code[j][1] * cosValue;
             
                // DFT_X_fd
            	factor = num_DFT*(i+1)*j;
            	cosValue = PETrigonometry.cos(factor);
				sinValue = PETrigonometry.sin(factor);	
				DFT_data[i+1][0] +=  X_fd[j][0] * cosValue + X_fd[j][1] * sinValue;
				DFT_data[i+1][1] += -X_fd[j][0] * sinValue + X_fd[j][1] * cosValue;
                // DFT_C
				DFT_code[i+1][0] +=  input_code[j][0] * cosValue + input_code[j][1] * sinValue;
                DFT_code[i+1][1] += -input_code[j][0] * sinValue + input_code[j][1] * cosValue;
                            
                // DFT_X_fd
            	factor = num_DFT*(i+2)*j;
            	cosValue = PETrigonometry.cos(factor);
				sinValue = PETrigonometry.sin(factor);	
				DFT_data[i+2][0] +=  X_fd[j][0] * cosValue + X_fd[j][1] * sinValue;
				DFT_data[i+2][1] += -X_fd[j][0] * sinValue + X_fd[j][1] * cosValue;
                // DFT_C
				DFT_code[i+2][0] +=  input_code[j][0] * cosValue + input_code[j][1] * sinValue;
                DFT_code[i+2][1] += -input_code[j][0] * sinValue + input_code[j][1] * cosValue;
                
                // DFT_X_fd
            	factor = num_DFT*(i+3)*j;
            	cosValue = PETrigonometry.cos(factor);
				sinValue = PETrigonometry.sin(factor);	
				DFT_data[i+3][0] +=  X_fd[j][0] * cosValue + X_fd[j][1] * sinValue;
				DFT_data[i+3][1] += -X_fd[j][0] * sinValue + X_fd[j][1] * cosValue;
                // DFT_C
				DFT_code[i+3][0] +=  input_code[j][0] * cosValue + input_code[j][1] * sinValue;
                DFT_code[i+3][1] += -input_code[j][0] * sinValue + input_code[j][1] * cosValue;      
                
            }
            /************************** step3: multiply_result calculation *************************************/
            multiply[i][0]= DFT_data[i][0]*DFT_code[i][0] - DFT_data[i][1]*(-DFT_code[i][1]);
            multiply[i][1]= DFT_code[i][0]*DFT_data[i][1] + DFT_data[i][0]*(-DFT_code[i][1]);
            
            multiply[i+1][0]= DFT_data[i+1][0]*DFT_code[i+1][0] - DFT_data[i+1][1]*(-DFT_code[i+1][1]);
            multiply[i+1][1]= DFT_code[i+1][0]*DFT_data[i+1][1] + DFT_data[i+1][0]*(-DFT_code[i+1][1]);
            
            multiply[i+2][0]= DFT_data[i+2][0]*DFT_code[i+2][0] - DFT_data[i+2][1]*(-DFT_code[i+2][1]);
            multiply[i+2][1]= DFT_code[i+2][0]*DFT_data[i+2][1] + DFT_data[i+2][0]*(-DFT_code[i+2][1]);
            
            multiply[i+3][0]= DFT_data[i+3][0]*DFT_code[i+3][0] - DFT_data[i+3][1]*(-DFT_code[i+3][1]);
            multiply[i+3][1]= DFT_code[i+3][0]*DFT_data[i+3][1] + DFT_data[i+3][0]*(-DFT_code[i+3][1]);
        }
		return multiply;
	}
	
	private void R_fd_Cal(float num_R_fd, float multiply[][]) {
		// TODO Auto-generated method stub
		float[][] R_fd = new float[nrOfSamples][2];
		float factor = 0.0f;
		float cosValue = 0.0f;
		float sinValue = 0.0f;
		float Rfd_square = 0.0f;
		for (int i=0; i<nrOfSamples; i=i+2) {
            for (int j=0; j<nrOfSamples; j=j+1) {
                factor = num_R_fd*i*j;
                cosValue = PETrigonometry.cos(factor);
                sinValue = PETrigonometry.sin(factor);
                R_fd[i][0] += multiply[j][0]*cosValue - multiply[j][1]*sinValue;
                R_fd[i][1] += multiply[j][0]*sinValue + multiply[j][1]*cosValue;
                
                factor = num_R_fd*(i+1)*j;
                cosValue = PETrigonometry.cos(factor);
                sinValue = PETrigonometry.sin(factor);
                R_fd[i+1][0] += multiply[j][0]*cosValue - multiply[j][1]*sinValue;
                R_fd[i+1][1] += multiply[j][0]*sinValue + multiply[j][1]*cosValue;
                
                factor = num_R_fd*(i+2)*j;
                cosValue = PETrigonometry.cos(factor);
                sinValue = PETrigonometry.sin(factor);
                R_fd[i+2][0] += multiply[j][0]*cosValue - multiply[j][1]*sinValue;
                R_fd[i+2][1] += multiply[j][0]*sinValue + multiply[j][1]*cosValue;
                
                factor = num_R_fd*(i+3)*j;
                cosValue = PETrigonometry.cos(factor);
                sinValue = PETrigonometry.sin(factor);
                R_fd[i+3][0] += multiply[j][0]*cosValue - multiply[j][1]*sinValue;
                R_fd[i+3][1] += multiply[j][0]*sinValue + multiply[j][1]*cosValue;
            }      
            R_fd[i][0] = R_fd[i][0]/nrOfSamples;
            R_fd[i][1] = R_fd[i][1]/nrOfSamples;
            Rfd_square = R_fd[i][0]*R_fd[i][0] + R_fd[i][1]*R_fd[i][1];
            if (Rfd_square >S_max) {
            	S_max = Rfd_square;
            	CodeVerschiebung = i;
			}
            R_fd[i+1][0] = R_fd[i+1][0]/nrOfSamples;
            R_fd[i+1][1] = R_fd[i+1][1]/nrOfSamples;
            Rfd_square = R_fd[i+1][0]*R_fd[i+1][0] + R_fd[i+1][1]*R_fd[i+1][1];
            if (Rfd_square >S_max) {
            	S_max = Rfd_square;
            	CodeVerschiebung = i+1;
			}
            R_fd[i+2][0] = R_fd[i+2][0]/nrOfSamples;
            R_fd[i+2][1] = R_fd[i+2][1]/nrOfSamples;
            Rfd_square = R_fd[i+2][0]*R_fd[i+2][0] + R_fd[i+2][1]*R_fd[i+2][1];
            if (Rfd_square >S_max) {
            	S_max = Rfd_square;
            	CodeVerschiebung = i+2;
			}
            R_fd[i+3][0] = R_fd[i+1][0]/nrOfSamples;
            R_fd[i+3][1] = R_fd[i+1][1]/nrOfSamples;
            Rfd_square = R_fd[i+3][0]*R_fd[i+3][0] + R_fd[i+3][1]*R_fd[i+3][1];
            if (Rfd_square >S_max) {
            	S_max = Rfd_square;
            	CodeVerschiebung = i+3;
			}
        }
	}

	public int getCodePhase() {
		// TODO Auto-generated method stub
		return CodeVerschiebung;
	}

}