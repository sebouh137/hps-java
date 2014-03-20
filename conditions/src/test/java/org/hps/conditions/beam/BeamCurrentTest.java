package org.hps.conditions.beam;

import static org.hps.conditions.ConditionsTableConstants.BEAM_CURRENT;

import java.io.File;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

import junit.framework.TestCase;

import org.hps.conditions.DatabaseConditionsManager;
import org.hps.conditions.beam.BeamCurrent.BeamCurrentCollection;
import org.lcsim.conditions.ConditionsManager;
import org.lcsim.event.EventHeader;
import org.lcsim.util.Driver;
import org.lcsim.util.cache.FileCache;
import org.lcsim.util.loop.LCSimLoop;

/**
 * This test checks the beam current values by run.
 * @author Jeremy McCormick <jeremym@slac.stanford.edu>
 */
public class BeamCurrentTest extends TestCase {
    
    /** This test file has a few events from the "good runs" of the Test Run. */
    private static final String TEST_FILE_URL = "http://www.lcsim.org/test/hps/conditions_test.slcio";
    
    /** Answer key for beam current by run. */
    static Map<Integer,Double> beamCurrentAnswerKey = new HashMap<Integer,Double>();
    
    /** Setup the beam current answer key by run. */
    static {
        beamCurrentAnswerKey.put(1349, 54879.7343788147); 
        beamCurrentAnswerKey.put(1351, 26928.0426635742);        
        beamCurrentAnswerKey.put(1353, 204325.132622242); 
        beamCurrentAnswerKey.put(1354, 148839.141475141); 
        beamCurrentAnswerKey.put(1358, 92523.9428218845); 
        beamCurrentAnswerKey.put(1359, 91761.4541434497);       
        beamCurrentAnswerKey.put(1360, 209883.979889035); 
        beamCurrentAnswerKey.put(1362, 110298.553449392); 
        beamCurrentAnswerKey.put(1363, 8556.8459701538); 
    }
    
    /**
     * Run the test.
     * @throws Exception 
     */
    public void test() throws Exception {

        // Cache file locally from URL.
        FileCache cache = new FileCache();
        File testFile = cache.getCachedFile(new URL(TEST_FILE_URL));
        
        // Create the LCSimLoop.
        LCSimLoop loop = new LCSimLoop();
        
        // Reconfigure the conditions system to override the manager created by LCSimLoop.
        DatabaseConditionsManager conditionsManager = DatabaseConditionsManager.createInstance();
        conditionsManager.configure("/org/hps/conditions/config/conditions_database_testrun_2013.xml");
        
        // Configure and run the loop.
        loop.setLCIORecordSource(testFile);
        loop.add(new BeamCurrentChecker());
        loop.loop(-1, null);
    }
    
    /**
     * This Driver will check the beam current for a run against the answer key.
     * @author Jeremy McCormick <jeremym@slac.stanford.edu>
     */
    class BeamCurrentChecker extends Driver {
        
        int currentRun = Integer.MIN_VALUE;
        
        /**
         * This method will check the beam current against the answer key
         * for the first event of a new run.
         */
        public void process(EventHeader event) {
            if (currentRun != event.getRunNumber()) {
                currentRun = event.getRunNumber();
                BeamCurrentCollection collection = ConditionsManager.defaultInstance().getCachedConditions(BeamCurrentCollection.class, BEAM_CURRENT).getCachedData();
                BeamCurrent beamCurrent = collection.get(0);
                System.out.println("Run " + event.getRunNumber() + " has integrated beam current " + beamCurrent.getIntegratedBeamCurrent() + " nC.");
                assertEquals("Wrong beam current for run.", beamCurrentAnswerKey.get(currentRun), beamCurrent.getIntegratedBeamCurrent());
            }
        }
    }
}
