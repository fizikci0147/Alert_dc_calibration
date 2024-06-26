/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import org.freehep.math.minuit.MnUserParameters;
import org.jlab.groot.math.Func1D;
import org.jlab.rec.dc.Constants;
/**
 *
 * @author ziegler
 */
public class FitLine extends Func1D{
    public int i;
    private FitFunction fc ;
    public FitLine() {
        super("fcn", 0.0, 2.0);
        fc = new FitFunction();
    }
    public static final int nPars = 11;
    private double[] par = new double[nPars];
    public FitLine(String name, int i, MnUserParameters pars) {
        super(name, 0.0, 2.0);
        this.i = i;
        fc = new FitFunction();
        this.initParameters(pars);
    }

    private void initParameters(MnUserParameters pars) {
        for(int p = 0; p< nPars; p++) {
            par[p] = pars.value(p);
        }
    }
    @Override
    public double evaluate(double x) { 
        double calcTime = 0;
        double B = 0;
        double ralpha=0;
        //local angle correction
        //reduce the corrected angle
        double v_0 = par[0];
        double vm = par[1];
        double tmax = par[3];
        double distbeta = par[4]; 
        double delBf = par[5]; 
        double Bb1 = par[6]; 
        double Bb2 = par[7]; 
        double Bb3 = par[8]; 
        double Bb4 = par[9]; 
        double R = par[2];
        double dmax = par[10];
        double deltatime_beta = (Math.sqrt(x * x + (distbeta * fc._beta * fc._beta) 
                * (distbeta* fc._beta * fc._beta)) - x) / Constants.V0AVERAGED;

        calcTime = fc.polyFcnMac(x,  ralpha,  B,  v_0,  vm,  R, 
            tmax,  dmax,  delBf,  Bb1,  Bb2,  Bb3,  Bb4, i+1) + deltatime_beta ;
        
        //System.out.println("ijk "+i+""+j+""+k+" b "+(float)T2DCalib.BfieldValues[k]+" ralpha "+(float)ralpha+" x "+x+" time "+(float)calcTime);
        return calcTime;
    }

    
}
