package org.hps.recon.tracking;

import java.io.File;
import java.net.URL;
import junit.framework.TestCase;
import org.hps.recon.tracking.gbl.GBLOutputDriver;
import org.hps.recon.tracking.gbl.GBLRefitterDriver;
import org.lcsim.util.cache.FileCache;
import org.lcsim.util.loop.LCSimLoop;

/**
 * This provides a template for testing track reconstruction issues
 *
 * @author Norman A Graf
 *
 * @version $Id:
 */
public class TrackRecoFromScratchTest extends TestCase
{
    static final String testURLBase = "http://www.lcsim.org/test/hps-java";
//    static final String testFileName = "TrackRecoFromScratchTest_radmuon_12.lcio-1-1788.slcio";
    static final String testFileName = "singleFullEnergyElectrons_SLIC-v05-00-00_Geant4-v10-01-02_QGSP_BERT_HPS-EngRun2015-Nominal-v2-fieldmap_1kEvents_recon_1Track_6Hits.slcio";
    private final int nEvents = 1;

    public void testRecon() throws Exception
    {
        File lcioInputFile = null;
        URL testURL = new URL(testURLBase + "/" + testFileName);
        FileCache cache = new FileCache();
        lcioInputFile = cache.getCachedFile(testURL);
        LCSimLoop loop = new LCSimLoop();
        loop.setLCIORecordSource(lcioInputFile);

        loop.add(new org.hps.recon.tracking.SimpleTrackerDigiDriver());
        loop.add(new org.hps.recon.tracking.HelicalTrackHitDriver());
        loop.add(new org.hps.recon.tracking.TrackerReconDriver());
        GBLRefitterDriver gblRefitter = new GBLRefitterDriver();
        loop.add(gblRefitter);
        GBLOutputDriver gbl = new org.hps.recon.tracking.gbl.GBLOutputDriver();
        gbl.setDebug(5);
        gbl.setIsMC(true);
        loop.add(gbl);

        try {
            loop.loop(nEvents);
        } catch (Exception e) {
            System.out.println("test should have failed");
            System.out.println("e");
        }
       
        loop.dispose();
    }
}