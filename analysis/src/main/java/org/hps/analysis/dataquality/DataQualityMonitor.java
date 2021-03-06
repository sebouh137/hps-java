package org.hps.analysis.dataquality;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.hps.conditions.beam.BeamEnergy.BeamEnergyCollection;
import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TIData;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * sort of an interface for DQM analysis drivers creates the DQM database
 * manager, checks whether row exists in db etc
 *
 * @author mgraham on Apr 15, 2014 update mgraham on May 15, 2014 to include
 * calculateEndOfRunQuantities & printDQMData i.e. useful methods
 */
public class DataQualityMonitor extends Driver {

    private static final Logger LOGGER = Logger.getLogger(DataQualityMonitor.class.getPackage().getName());

    protected Double beamEnergy;
    public void setBeamEnergy(double e){
    this.beamEnergy = e;
    }
    public double getBeamEnergy(){
    return this.beamEnergy;
    }
    
    protected AIDA aida = AIDA.defaultInstance();
    protected DQMDatabaseManager manager;
    protected String recoVersion = "v0.0";
    protected static int runNumber = 1350;
    protected boolean overwriteDB = false;
    protected boolean connectToDB = false;
    protected boolean printDQMStrings = false;
    protected Map<String, Double> monitoredQuantityMap = new HashMap<>();
    protected boolean debug = false;
    protected boolean outputPlots = false;
    protected String outputPlotDir = "DQMOutputPlots/";

    @Override
    protected void detectorChanged(Detector detector){
        BeamEnergyCollection beamEnergyCollection = 
            this.getConditionsManager().getCachedConditions(BeamEnergyCollection.class, "beam_energies").getCachedData();        
        if(beamEnergy== null && beamEnergyCollection != null && beamEnergyCollection.size() != 0)
            beamEnergy = beamEnergyCollection.get(0).getBeamEnergy();
        else{
            LOGGER.log(Level.WARNING, "warning:  beam energy not found.  Using a 6.6 GeV as the default energy");
            beamEnergy = 6.6;
        }
       
    }
    
    
    String triggerType = "all";//allowed types are "" (blank) or "all", singles0, singles1, pairs0,pairs1
    public boolean isGBL = false;

    public void setTriggerType(String type) {
        this.triggerType = type;
    }

    public void setIsGBL(boolean isgbl) {
        this.isGBL = isgbl;
    }

    public void setRecoVersion(String recoVersion) {
        this.recoVersion = recoVersion;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public void setRunNumber(int run) {
        DataQualityMonitor.runNumber = run;
    }

    public void setOverwriteDB(boolean overwrite) {
        this.overwriteDB = overwrite;
    }

    public void setConnectToDB(boolean connect) {
        this.connectToDB = connect;
    }

    public void setPrintDQMStrings(boolean print) {
        this.printDQMStrings = print;
    }

    public void setOutputPlots(boolean out) {
        this.outputPlots = out;
    }

    public void setOutputPlotDir(String dir) {
        this.outputPlotDir = dir;
    }

    @Override
    public void endOfData() {
        calculateEndOfRunQuantities();
        fillEndOfRunPlots();
        printDQMData();
        if (printDQMStrings) {
            printDQMStrings();
        }
        LOGGER.info("Write to database =  " + connectToDB);
        if (connectToDB) {
            LOGGER.info("Connecting To Database...getting DQMDBManager");
            manager = DQMDatabaseManager.getInstance();
            //check to see if I need to make a new db entry
            boolean entryExists = false;
            try {
                entryExists = checkRowExists();
                if (entryExists) {
                    LOGGER.info("Found an existing run/reco entry in the dqm database; overwrite = " + overwriteDB);
                }
            } catch (SQLException ex) {
                Logger.getLogger(DataQualityMonitor.class.getName()).log(Level.SEVERE, null, ex);
            }

            if (!entryExists) {
                makeNewRow();
            }
            dumpDQMData();
        }
    }

    private void makeNewRow() {
        LOGGER.info("is the data base connected?  " + manager.isConnected);
        if (manager.isConnected) {
            String ins = "insert into dqm SET runNumber=" + runNumber;
//            LOGGER.info(ins);
            manager.updateQuery(ins);
            ins = "update  dqm SET recoVersion='" + recoVersion + "' where runNumber=" + runNumber;
            manager.updateQuery(ins);
        }
        LOGGER.info("Made a new row for runNumber=" + runNumber + "; recoVersion=" + recoVersion);
    }

    private boolean checkRowExists() throws SQLException {
        String ins = "select * from dqm where " + getRunRecoString();
        ResultSet res = manager.selectQuery(ins);
        return res.next(); //this is a funny way of determining if the ResultSet has any entries
    }

    public boolean checkSelectionIsNULL(String var) throws SQLException {
        String ins = "select " + var + " from dqm where " + getRunRecoString();
        ResultSet res = manager.selectQuery(ins);
        res.next();
        double result = res.getDouble(var);
        if (res.wasNull()) {
            return true;
        }
        LOGGER.info("checkSelectionIsNULL::" + var + " = " + result);
        return false;
    }

    public String getRunRecoString() {
        return "runNumber=" + runNumber + " and recoVersion='" + recoVersion + "'";
    }

    //override this method to do something interesting   
    //like fill some plots that you only want to fill at end of data (e.g. for occupancies)
    public void fillEndOfRunPlots() {
    }

    //override this method to do something interesting   
    //like calculate averages etc. that can then be put in the db  
    public void calculateEndOfRunQuantities() {
    }

    public void dumpDQMData() {
        for (Map.Entry<String, Double> entry : monitoredQuantityMap.entrySet()) {
            String name = entry.getKey();
            double val = entry.getValue();
            boolean isnull = false;
            try {
                isnull = checkSelectionIsNULL(name);
            } catch (SQLException ex) {
                Logger.getLogger(SvtMonitoring.class.getName()).log(Level.SEVERE, null, ex);
            }
            if (!overwriteDB && !isnull) {
                LOGGER.info("Not writing because " + name + " is already filled for this entry");
                continue; //entry exists and I don't want to overwrite                
            }
            String put = "update dqm SET " + name + " = " + val + " WHERE " + getRunRecoString();
            LOGGER.info(put);
            manager.updateQuery(put);

        }
    }

    public String getTriggerType() {
        return triggerType;
    }

    public boolean matchTriggerType(TIData triggerData) {
        if (triggerType.contentEquals("") || triggerType.contentEquals("all")) {
            return true;
        }
        if (triggerData.isSingle0Trigger() && triggerType.contentEquals("singles0")) {
            return true;
        }
        if (triggerData.isSingle1Trigger() && triggerType.contentEquals("singles1")) {
            return true;
        }
        if (triggerData.isPair0Trigger() && triggerType.contentEquals("pairs0")) {
            return true;
        }
        if (triggerData.isPair1Trigger() && triggerType.contentEquals("pairs1")) {
            return true;
        }
        return false;

    }

    public boolean matchTrigger(EventHeader event) {
        boolean match = true;
        if (event.hasCollection(GenericObject.class, "TriggerBank")) {
            List<GenericObject> triggerList = event.get(GenericObject.class, "TriggerBank");
            for (GenericObject data : triggerList) {
                if (AbstractIntData.getTag(data) == TIData.BANK_TAG) {
                    TIData triggerData = new TIData(data);
                    if (!matchTriggerType(triggerData))//only process singles0 triggers...
                    {
                        match = false;
                    }
                }
            }
        } else if (debug) {
            LOGGER.info(this.getClass().getSimpleName() + ":  No trigger bank found...running over all trigger types");
        }
        return match;
    }

    //override this method to do something interesting   
    //like print the DQM data log file
    public void printDQMData() {
    }

    //override this method to do something interesting   
    //like print the DQM db variable strings in a good 
    //format for making the db column headers
    public void printDQMStrings() {
    }
}
