package org.clas.detector.clas12calibration.dc.t2d;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.timetodistance.T2DFunctions;
import org.jlab.utils.groups.IndexedTable;


public class TableLoader {

    public TableLoader() {
            // TODO Auto-generated constructor stub
    }
    static final protected int nBinsT=2000;
    //public static double[][][][][] DISTFROMTIME = new double[6][6][6][6][850]; // sector slyr alpha Bfield time bins
    
    static boolean T2DLOADED = false;
    static boolean T0LOADED = false;
    
    static int minBinIdxT  = 0;
    static int[][][][] maxBinIdxT  = new int[6][6][8][6];
    public static double[][][] DISTFROMTIME = new double[6][6][nBinsT]; // sector slyr time bins [s][r][ibfield][icosalpha][tbin]
    
    
    
    public static synchronized void FillT0Tables(int run, String variation) {
        if (T0LOADED) return;
        System.out.println(" T0 TABLE FILLED..... for Run "+run+" with VARIATION "+variation);
        DatabaseConstantProvider dbprovider = new DatabaseConstantProvider(run, variation);
        dbprovider.loadTable("/calibration/dc/time_corrections/T0Corrections");
        //disconnect from database. Important to do this after loading tables.
        dbprovider.disconnect();
        // T0-subtraction
        double[][][][] T0 ;
        double[][][][] T0ERR ;
        //T0s
        T0 = new double[6][6][7][6]; //nSec*nSL*nSlots*nCables
        T0ERR = new double[6][6][7][6]; //nSec*nSL*nSlots*nCables
        for (int i = 0; i < dbprovider.length("/calibration/dc/time_corrections/T0Corrections/Sector"); i++) {
            int iSec = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Sector", i);
            int iSly = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Superlayer", i);
            int iSlot = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Slot", i);
            int iCab = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Cable", i);
            double t0 = dbprovider.getDouble("/calibration/dc/time_corrections/T0Corrections/T0Correction", i);
            double t0Error = dbprovider.getDouble("/calibration/dc/time_corrections/T0Corrections/T0Error", i);

            T0[iSec - 1][iSly - 1][iSlot - 1][iCab - 1] = t0; 
            T0ERR[iSec - 1][iSly - 1][iSlot - 1][iCab - 1] = t0Error;
            Constants.setT0(T0);
            Constants.setT0Err(T0ERR);
            //System.out.println("T0 = "+t0);
        }
        T0LOADED = true;
    }
 
    public static int maxTBin = -1;

    public static synchronized void Fill(IndexedTable tab) {
        //CCDBTables 0 =  "/calibration/dc/signal_generation/doca_resolution";
        //CCDBTables 1 =  "/calibration/dc/time_to_distance/t2d";
        //CCDBTables 2 =  "/calibration/dc/time_corrections/T0_correction";	
        if (T2DLOADED) return;
        
        double stepSize = 0.0010;
        DecimalFormat df = new DecimalFormat("#");
        df.setRoundingMode(RoundingMode.CEILING);
        
        for(int s = 0; s<6; s++ ){ // loop over sectors

                for(int r = 0; r<6; r++ ){ //loop over slys
                    // Fill constants
                    delta_T0[s][r] = tab.getDoubleValue("delta_T0", s+1,r+1,0);
                    FracDmaxAtMinVel[s][r] = tab.getDoubleValue("c1", s+1,r+1,0);//use same table. names strings 
                    deltanm[s][r] = tab.getDoubleValue("deltanm", s+1,r+1,0);
                    v0[s][r] = tab.getDoubleValue("v0", s+1,r+1,0);
                    vmid[s][r] = tab.getDoubleValue("c2", s+1,r+1,0);
                    delta_bfield_coefficient[s][r] = tab.getDoubleValue("delta_bfield_coefficient", s+1,r+1,0); 
                    distbeta[s][r] = tab.getDoubleValue("distbeta", s+1,r+1,0); 
                    b1[s][r] = tab.getDoubleValue("b1", s+1,r+1,0);
                    b2[s][r] = tab.getDoubleValue("b2", s+1,r+1,0);
                    b3[s][r] = tab.getDoubleValue("b3", s+1,r+1,0);
                    b4[s][r] = tab.getDoubleValue("b4", s+1,r+1,0);
                    Tmax[s][r] = tab.getDoubleValue("tmax", s+1,r+1,0);
                    // end fill constants
                    double dmax = 2.*Constants.wpdist[r]; 

                                    
                }
        }	
        TableLoader.fillMissingTableBins();
        //TableLoader.test();
        System.out.println(" T2D TABLE FILLED.....");
        T2DLOADED = true;
     }
    
    
    public static synchronized void ReFill() {
        //reset
        DISTFROMTIME = new double[6][6][nBinsT]; // sector slyr alpha Bfield time bins [s][r][ibfield][icosalpha][tbin]
        minBinIdxT  = 0;
        maxBinIdxT  = new int[6][6][8][6];
        maxTBin = -1;
        double stepSize = 0.0010;
        DecimalFormat df = new DecimalFormat("#");
        df.setRoundingMode(RoundingMode.CEILING);
        

        TableLoader.fillMissingTableBins();
        //TableLoader.test();
        System.out.println(" T2D TABLE RE-FILLED.....");
     }
    

    private static void fillMissingTableBins() {
        
        for(int s = 0; s<6; s++ ){ // loop over sectors

            for(int r = 0; r<6; r++ ){ //loop over slys
                    
                        
                        for(int tbin = 0; tbin<maxTBin; tbin++) {
                            if(DISTFROMTIME[s][r][tbin]!=0 && DISTFROMTIME[s][r][tbin+1]==0) {
                                DISTFROMTIME[s][r][tbin+1] = DISTFROMTIME[s][r][tbin];
                            }
                        }
                        
            }
        }
    }
    /**
     * 
     * @param x distance to wire in cm
     * @param alpha local angle in deg
     * @param bfield B field value a x in T
     * @param sector sector  
     * @param superlayer superlayer 
     * @return returns time (ns) when given inputs of distance x (cm), local angle alpha (degrees) and magnitude of bfield (Tesla).  
     */
    public static synchronized double calc_Time(double x, int sector, int superlayer) {
        int s = sector - 1;
        int r = superlayer - 1;
        double dmax = 2.*Constants.wpdist[r]; 
        double tmax = Tmax[s][r];
        double Bb1 = b1[s][r];
        double Bb2 = b2[s][r];
        double Bb3 = b3[s][r];
        double Bb4 = b4[s][r];
        if(x>dmax)
            x=dmax;
        double alpha=25; //we dont need these two, just keeping it for now so that code compiles
        double delBf =0.5;
        double bfield=0.5;
        return T2DFunctions.polyFcnMac(x,alpha,bfield, v0[s][r], vmid[s][r], FracDmaxAtMinVel[s][r], 
                tmax, dmax,delBf, Bb1, Bb2, Bb3, Bb4, superlayer) ;
        
    }
    
    public static double[][] delta_T0 = new double[6][6];
    public static double[][] delta_bfield_coefficient = new double[6][6];
    public static double[][] distbeta = new double[6][6];
    public static double[][] deltanm = new double[6][6];
    public static double[][] vmid = new double[6][6];
    public static double[][] v0 = new double[6][6];
    public static double[][] b1 = new double[6][6];
    public static double[][] b2 = new double[6][6];
    public static double[][] b3 = new double[6][6];
    public static double[][] b4 = new double[6][6];
    public static double[][] Tmax = new double[6][6];
    public static double[][] FracDmaxAtMinVel = new double[6][6];		// fraction of dmax corresponding to the point in the cell where the velocity is minimal

}
