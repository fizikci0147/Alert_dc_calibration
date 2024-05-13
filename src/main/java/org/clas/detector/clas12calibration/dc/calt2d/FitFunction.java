/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.freehep.math.minuit.FCNBase;
import org.jlab.groot.data.GraphErrors;
import org.jlab.rec.dc.Constants;

/**
 *
 * @author ziegler
 */
public class FitFunction implements FCNBase{

    public double _beta = 1.0;
    
    private Map<Coordinate, GraphErrors> _tvstrkdocasProf;
    private int i;
    
    public FitFunction() {
        
    }
    public FitFunction(int i, Map<Coordinate, GraphErrors> tvstrkdocasProf) {
        this.i = i;
        _tvstrkdocasProf = tvstrkdocasProf;
    }
         
    public double eval(double x, double[] par) {
        double ralpha=0., B=0.;
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
        double deltatime_beta = (Math.sqrt(x * x + (distbeta * _beta * _beta) 
                * (distbeta* _beta * _beta)) - x) / Constants.V0AVERAGED;

        double calcTime = this.polyFcnMac(x,  ralpha,  B,  v_0,  vm,  R, 
            tmax,  dmax,  delBf,  Bb1,  Bb2,  Bb3,  Bb4, i+1) + deltatime_beta ;
        
        return calcTime;
    }
    public double polyFcnMac(double x, double alpha, double bfield, double v_0, double vm, double R, 
            double tmax, double dmax, double delBf, double Bb1, double Bb2, double Bb3, double Bb4, int superlayer) {
        
        if(x>dmax)
            x=dmax;
        double time = 0;
        //   rcapital is an intermediate parameter
        double rcapital = R*dmax;
        //   delt is another intermediate parameter
        double delt=tmax-dmax/v_0;
        double delv=1./vm-1./v_0;
        //   now calculate the primary parameters a, b, c, d
        
        double a=0.1,b=0.1,c=0.1,d=0.1;
        time = a*x*x*x*x + b*x*x*x + c*x*x + d*x ;
        
        return time;
    }
    @Override
    public double valueOf(double[] par) {
        double chisq = 0;
        double delta = 0;
                    if(_tvstrkdocasProf.get(new Coordinate(this.i)).getVectorX().size()>0){ 
                        GraphErrors gr = _tvstrkdocasProf.get(new Coordinate(this.i));
                            
                        for (int ix =0; ix< gr.getDataSize(0); ix++) {
                            double x = gr.getDataX(ix);
                            double time = gr.getDataY(ix);
                            double err = gr.getDataEY(ix);
                            if(err>0) {
                                double calcTime = this.eval(x, par);
                                delta = (time - calcTime) / err; 
                                chisq += delta * delta;
                            }
                        }
                    }
        return chisq;
        
    }
    
    
}
