package org.hps.analysis.dataquality;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * sort of an interface for DQM analysis drivers
 * creates the DQM database manager, checks whether row exists in db etc
 * @author mgraham on Apr 15, 2014
 */
public class DataQualityMonitor extends Driver {

    public AIDA aida = AIDA.defaultInstance();
    public DQMDatabaseManager manager;
    public String recoVersion = "v0.0";
    public static int runNumber = 1350;
    public boolean overwriteDB = false;

    public void setRecoVersion(String recoVersion) {
        this.recoVersion = recoVersion;
    }

    public void setOverwriteDB(boolean overwrite) {
        this.overwriteDB = overwrite;
    }

    public void DataQualityMonitor() {

    }

    public void endOfData() {
        manager = DQMDatabaseManager.getInstance();
        //fill any plots that only get filled at end of data...e.g. SVT occupancy plots
        fillEndOfRunPlots();

        //check to see if I need to make a new db entry
        boolean entryExists = false;
        try {
            entryExists = checkRowExists();
            if(entryExists)System.out.println("Found an existing run/reco entry in the dqm database");
        } catch (SQLException ex) {
            Logger.getLogger(DataQualityMonitor.class.getName()).log(Level.SEVERE, null, ex);
        }

        if (!entryExists)
            makeNewRow();
        else
            if (!overwriteDB)
                return; //entry exists and I don't want to overwrite
        
        dumpDQMData();
    }

    private void makeNewRow() {
        System.out.println("is the data base connected?  " + manager.isConnected);
        if (manager.isConnected) {
            String ins = "insert into dqm SET run=" + runNumber;
//            System.out.println(ins);
            manager.updateQuery(ins);
            ins = "update  dqm SET recoversion='" + recoVersion + "' where run=" + runNumber;
            manager.updateQuery(ins);
        }
        System.out.println("Made a new row for run=" + runNumber + "; recon version=" + recoVersion);
    }

    private boolean checkRowExists() throws SQLException {
        String ins = "select * from dqm where " + getRunRecoString();
        ResultSet res = manager.selectQuery(ins);
        if (res.next()) //this is a funny way of determining if the ResultSet has any entries
            return true;
        return false;
    }

    public String getRunRecoString() {
        return "run=" + runNumber + " and recoversion='"+recoVersion+"'";
    }

    //override this method to do something interesting   
    public void fillEndOfRunPlots() {
    }
    //override this method to do something interesting   
    public void dumpDQMData() {
    }

}
