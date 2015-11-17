package org.hps.users.spaul.feecc;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class CustomBinning {
	public CustomBinning(File f) throws FileNotFoundException{
		Scanner s = new Scanner(f);

		nTheta = s.nextInt();  //number of bins in theta;
		thetaMax = new double[nTheta];
		thetaMin = new double[nTheta];
		
		phiMax = new double[nTheta][];
		phiMin = new double[nTheta][];
		int i = 0;
		while(s.hasNext()){ //new row
			int nPhi = s.nextInt();
			thetaMin[i] = s.nextDouble();
			thetaMax[i] = s.nextDouble();
			
			phiMax[i] = new double[nPhi];
			phiMin[i] = new double[nPhi];
			for(int j = 0; j<nPhi; j++){
				phiMin[i][j] = s.nextDouble();
				phiMax[i][j] = s.nextDouble();
			}
			i++;
		}
	}
	double[][] phiMax;
	double[][] phiMin;
	public double thetaMax[], thetaMin[];
	public int nTheta;

	double getSteradians(int binNumber){
		double t1 = thetaMin[binNumber];
		double t2 = thetaMax[binNumber];
		double dCos = Math.cos(t1)-Math.cos(t2);
		double dPhiTot = 0;
		for(int i = 0; i< phiMax[binNumber].length; i++){
			dPhiTot += phiMax[binNumber][i]-phiMin[binNumber][i];
		}
		return 2*dPhiTot*dCos;  //factor of two because top and bottom
	}
	boolean inRange(double theta, double phi){
		phi = Math.abs(phi);
		/*int i =(int) Math.floor((theta-theta0)/deltaTheta);
		if(i>= nTheta || i<0)
			return false;*/
		if(theta > thetaMax[nTheta-1] || theta < thetaMin[0])
			return false;
		int i;
		boolean found = false;
		for(i = 0; i< nTheta; i++){
			if(theta > thetaMin[i] && theta < thetaMax[i]){
				found = true;
				break;
			}
		}
		if(!found)
			return false;
		
		for(int j = 0; j<phiMax[i].length; j++){
			if(phi>phiMin[i][j] && phi< phiMax[i][j])
				return true;
		}
		return false;

	}
	public double getTotSteradians() {
		double tot = 0;
		for(int i = 0; i<nTheta; i++){
			tot += getSteradians(i);
		}
		return tot;
	}
	/**
	 * @param bin
	 * @param a = 2E/M
	 * @return the integral of 1/sin^4(th/2)*cos^2(th/2)/(1+a*sin^2(th/2)) times dPhi,
	 * which appears in the integral of mott scattering.
	 */
	public double mottIntegralFactor(double a, int bin){
		double dPhi = 0;
		for(int i = 0; i< phiMax[bin].length; i++)
			dPhi += 2*(phiMax[bin][0] - phiMin[bin][0]); //factor of 2 from top and bottom
		
		double csc2 = Math.pow(Math.sin(thetaMax[bin]/2), -2);
		double Imax = (-csc2+(1+a)*Math.log(a+2*csc2));
		       csc2 = Math.pow(Math.sin(thetaMin[bin]/2), -2);
		double Imin = (-csc2+(1+a)*Math.log(a+2*csc2));
		return 2*dPhi*(Imax-Imin);
	}
}
