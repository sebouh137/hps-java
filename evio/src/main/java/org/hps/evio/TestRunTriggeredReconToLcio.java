package org.hps.evio;

import hep.physics.event.generator.MCEvent;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.readout.ecal.ClockSingleton;
import org.hps.readout.ecal.ReadoutTimestamp;
import org.hps.readout.ecal.TriggerableDriver;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.base.BaseLCSimEvent;
import org.lcsim.geometry.Detector;
import org.lcsim.lcio.LCIOWriter;

/**
 * This class takes raw data generated from MC and converts it to EVIO. The goal
 * is to make this look like data which will come off the actual ET ring during
 * the test run.
 */
public class TestRunTriggeredReconToLcio extends TriggerableDriver {

    String rawCalorimeterHitCollectionName = "EcalReadoutHits";
    String outputFile = "TestRunData.slcio";
    private int eventsWritten = 0;
    private int eventNum = 0;
    //interval for trigger candidates (tridents, A'), if used
    private int triggerSpacing = 250;
    private boolean rejectBackground = false;
//    HPSEcalConditions ecalIDConverter = null;
    EcalHitWriter ecalWriter = null;
    SVTHitWriter svtWriter = null;
    TriggerDataWriter triggerWriter = null;
    List<HitWriter> writers = null;
    LCIOWriter lcioWriter = null;
    Queue<EventHeader> events = null;
    private int ecalMode = EventConstants.ECAL_PULSE_INTEGRAL_MODE;
    List<MCParticle> mcParticles = null;
    List<SimTrackerHit> trackerHits = null;
    List<SimCalorimeterHit> ecalHits = null;
    List<SimTrackerHit> ecalScoringPlaneHits = null;
    //MC collections from the last 500n'th event (trident or preselected trigger event)
    List<MCParticle> triggerMCParticles = null;
    List<SimTrackerHit> triggerTrackerHits = null;
    List<SimCalorimeterHit> triggerECalHits = null;
    List<SimTrackerHit> triggerECalScoringPlaneHits = null;
    static final String ecalCollectionName = "EcalHits";
    static final String trackerCollectionName = "TrackerHits";
    private final String relationCollectionName = "SVTTrueHitRelations";
    String ecalScoringPlaneHitsCollectionName = "TrackerHitsECal";
    private int verbosity = 1;
    private boolean writeSvtData = true;
    private boolean writeEcalData = true;
    private boolean writeTriggerData = true;

    public TestRunTriggeredReconToLcio() {
        setTriggerDelay(0);
    }

    @Override
    public void detectorChanged(Detector detector) {
        // set the detector
        ecalWriter.setDetector(detector);
    }

    public void setEcalMode(int ecalMode) {
        this.ecalMode = ecalMode;
        if (ecalMode != EventConstants.ECAL_RAW_MODE && ecalMode != EventConstants.ECAL_PULSE_MODE && ecalMode != EventConstants.ECAL_PULSE_INTEGRAL_MODE) {
            throw new IllegalArgumentException("invalid mode " + ecalMode);
        }
        if (ecalWriter != null) {
            ecalWriter.setMode(ecalMode);
        }
    }

    public void setOutputFile(String outputFile) {
        this.outputFile = outputFile;
    }

    public void setTriggerSpacing(int triggerSpacing) {
        this.triggerSpacing = triggerSpacing;
    }

    public void setRejectBackground(boolean rejectBackground) {
        this.rejectBackground = rejectBackground;
    }

    public void setRawCalorimeterHitCollectionName(String rawCalorimeterHitCollectionName) {
        this.rawCalorimeterHitCollectionName = rawCalorimeterHitCollectionName;
        if (ecalWriter != null) {
            ecalWriter.setHitCollectionName(rawCalorimeterHitCollectionName);
        }
    }

    /**
     * Set the amount of printouts generated by the writers. 0 = silent, 1 =
     * normal, 2+ = debug
     *
     * @param verbosity
     */
    public void setVerbosity(int verbosity) {
        this.verbosity = verbosity;
        if (writers != null) {
            for (HitWriter hitWriter : writers) {
                hitWriter.setVerbosity(verbosity);
            }
        }
    }

    /**
     * Set whether the LCIO writer looks for SVT readout data.
     *
     * @param writeSvtData True by default.
     */
    public void setWriteSvtData(boolean writeSvtData) {
        this.writeSvtData = writeSvtData;
    }

    /**
     * Set whether the LCIO writer looks for ECal readout data.
     *
     * @param writeEcalData True by default.
     */
    public void setWriteEcalData(boolean writeEcalData) {
        this.writeEcalData = writeEcalData;
    }

    /**
     * Set whether the LCIO writer looks for trigger readout data.
     *
     * @param writeTriggerData True by default.
     */
    public void setWriteTriggerData(boolean writeTriggerData) {
        this.writeTriggerData = writeTriggerData;
    }

    @Override
    protected void startOfData() {
        super.startOfData();
        writers = new ArrayList<HitWriter>();

        if (writeEcalData) {
            ecalWriter = new EcalHitWriter();
            ecalWriter.setMode(ecalMode);
            ecalWriter.setHitCollectionName(rawCalorimeterHitCollectionName);
            writers.add(ecalWriter);
        }

        if (writeSvtData) {
            svtWriter = new SVTHitWriter();
            writers.add(svtWriter);
        }

        if (writeTriggerData) {
            triggerWriter = new TriggerDataWriter();
            writers.add(triggerWriter);
        }

        for (HitWriter hitWriter : writers) {
            hitWriter.setVerbosity(verbosity);
        }

        try {
            lcioWriter = new LCIOWriter(new File(outputFile));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        events = new LinkedList<EventHeader>();
    }

    @Override
    protected void endOfData() {
        if (verbosity >= 1) {
            System.out.println(this.getClass().getSimpleName() + " - wrote " + eventsWritten + " events in job; " + events.size() + " incomplete events in queue.");
        }
        try {
            lcioWriter.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    protected void process(EventHeader event) {
        if (event.hasCollection(SimCalorimeterHit.class, ecalCollectionName) && !event.get(SimCalorimeterHit.class, ecalCollectionName).isEmpty()) {
            mcParticles = event.getMCParticles();
            ecalHits = event.get(SimCalorimeterHit.class, ecalCollectionName);
            if (event.hasCollection(SimTrackerHit.class, trackerCollectionName)) {
                trackerHits = event.get(SimTrackerHit.class, trackerCollectionName);
            }
            if (event.hasCollection(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName)) {
                ecalScoringPlaneHits = event.get(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName);
            }
        }
        if (ClockSingleton.getClock() % triggerSpacing == 0) {
            if (event.hasCollection(MCParticle.class)) {
                triggerMCParticles = event.getMCParticles();
                triggerECalHits = event.get(SimCalorimeterHit.class, ecalCollectionName);
                if (event.hasCollection(SimTrackerHit.class, trackerCollectionName)) {
                    triggerTrackerHits = event.get(SimTrackerHit.class, trackerCollectionName);
                }
                if (event.hasCollection(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName)) {
                    triggerECalScoringPlaneHits = event.get(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName);
                }
            } else {
                triggerMCParticles = null;
                triggerECalHits = null;
                triggerTrackerHits = null;
                triggerECalScoringPlaneHits = null;
            }
        }

        checkTrigger(event);

        writerLoop:
        for (HitWriter hitWriter : writers) {
            if (hitWriter.hasData(event)) {
                if (verbosity >= 1) {
                    System.out.println(hitWriter.getClass().getSimpleName() + ": writing data, event " + event.getEventNumber());
                }
                for (EventHeader queuedEvent : events) {
                    if (!hitWriter.hasData(queuedEvent)) {
                        // Write data.
                        hitWriter.writeData(event, queuedEvent);
                        continue writerLoop;
                    }
                }

                throw new RuntimeException("no queued events waiting for an " + hitWriter.getClass().getSimpleName() + " bank");
            }
        }

        eventLoop:
        while (!events.isEmpty()) {
            EventHeader queuedEvent = events.peek();
            for (HitWriter hitWriter : writers) {
                if (!hitWriter.hasData(queuedEvent)) {
                    break eventLoop;
                }
            }
            events.poll();

            boolean writeThisEvent = true;
            if (rejectBackground && queuedEvent.hasCollection(LCRelation.class, relationCollectionName)) {
                writeThisEvent = false;
                List<LCRelation> trueHitRelations = queuedEvent.get(LCRelation.class, relationCollectionName);
                List<SimTrackerHit> trueHits = queuedEvent.get(SimTrackerHit.class, trackerCollectionName);
                for (LCRelation relation : trueHitRelations) {
                    if (trueHits.contains((SimTrackerHit) relation.getTo())) {
                        writeThisEvent = true;
                        break;
                    }
                }
            }

            if (writeThisEvent) {
                // Write this event.
                if (verbosity >= 1) {
                    System.out.println("writing filled LCIO event, event " + queuedEvent.getEventNumber());
                }
                try {
                    lcioWriter.write(queuedEvent);
                    ++eventsWritten;
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            } else {
                if (verbosity >= 1) {
                    System.out.println("rejecting filled LCIO event, event " + queuedEvent.getEventNumber() + " contains no SVT hits from truth particles");
                }
            }
        }
    }

    @Override
    protected void processTrigger(EventHeader event) {
        // Create an LCSim event and pass a flag so that conditions updates are disabled. --JM
        EventHeader lcsimEvent = new BaseLCSimEvent(DatabaseConditionsManager.getInstance().getRun(), event.getEventNumber(), event.getDetectorName(), (long) 4 * (Math.round(ClockSingleton.getTime() / 4)), false);
        events.add(lcsimEvent);
        if (verbosity >= 1) {
            System.out.println("Creating LCIO event " + eventNum);
        }
        if (triggerMCParticles == null || triggerMCParticles.isEmpty()) {
            lcsimEvent.put(MCEvent.MC_PARTICLES, mcParticles);
            if (verbosity >= 1) {
                System.out.println("Adding " + mcParticles.size() + " MCParticles");
            }
            if (ecalHits != null) {
                lcsimEvent.put(ecalCollectionName, ecalHits, SimCalorimeterHit.class, 0xe0000000);
                if (verbosity >= 1) {
                    System.out.println("Adding " + ecalHits.size() + " SimCalorimeterHits");
                }
            }
            if (trackerHits != null) {
                lcsimEvent.put(trackerCollectionName, trackerHits, SimTrackerHit.class, 0xc0000000);
                if (verbosity >= 1) {
                    System.out.println("Adding " + trackerHits.size() + " SimTrackerHits");
                }
            }
            if (ecalScoringPlaneHits != null) {
                lcsimEvent.put(ecalScoringPlaneHitsCollectionName, ecalScoringPlaneHits, SimTrackerHit.class, 0xc0000000);
                if (verbosity >= 1) {
                    System.out.println("Adding " + ecalScoringPlaneHits.size() + " ECalTrackerHits");
                }
            }
        } else {
            lcsimEvent.put(MCEvent.MC_PARTICLES, triggerMCParticles);
            if (verbosity >= 1) {
                System.out.println("Adding " + triggerMCParticles.size() + " MCParticles");
            }
            if (triggerECalHits != null) {
                lcsimEvent.put(ecalCollectionName, triggerECalHits, SimCalorimeterHit.class, 0xe0000000);
                if (verbosity >= 1) {
                    System.out.println("Adding " + triggerECalHits.size() + " SimCalorimeterHits");
                }
            }
            if (triggerTrackerHits != null) {
                lcsimEvent.put(trackerCollectionName, triggerTrackerHits, SimTrackerHit.class, 0xc0000000);
                if (verbosity >= 1) {
                    System.out.println("Adding " + triggerTrackerHits.size() + " SimTrackerHits");
                }
            }
            if (triggerECalScoringPlaneHits != null) {
                lcsimEvent.put(ecalScoringPlaneHitsCollectionName, triggerECalScoringPlaneHits, SimTrackerHit.class, 0xc0000000);
                if (verbosity >= 1) {
                    System.out.println("Adding " + triggerECalScoringPlaneHits.size() + " ECalTrackerHits");
                }
            }
        }
        lcsimEvent.put(ReadoutTimestamp.collectionName, event.get(ReadoutTimestamp.class, ReadoutTimestamp.collectionName));
        ++eventNum;
    }

    @Override
    public int getTimestampType() {
        return ReadoutTimestamp.SYSTEM_TRIGGERTIME;
    }
}
