package org.hps.analysis.dataquality;

import hep.aida.IHistogramFactory;
import hep.aida.IProfile1D;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.hps.analysis.examples.TrackAnalysis;
import org.hps.recon.tracking.FindableTrack;
import org.hps.recon.tracking.FittedRawTrackerHit;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.MCParticle;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.fit.helicaltrack.HelicalTrackCross;
import org.lcsim.fit.helicaltrack.HelixParamCalculator;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.IDDecoder;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHit;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;

/**
 *  DQM driver for the monte carlo track efficiency; makes a bunch of efficiency vs variable plots
 *  for all tracks and just electrons from trident/A' event, as well as "findable" tracks
 *  use the debugTrackEfficiency flag to print out info regarding individual failed events
 * @author mgraham on Mar 28, 2014
 */
// TODO:  Add some quantities for DQM monitoring:  e.g. <efficiency>, <eff>_findable
public class TrackMCEfficiency extends DataQualityMonitor {

    private String rawTrackerHitCollectionName = "SVTRawTrackerHits";
    private String helicalTrackHitCollectionName = "HelicalTrackHits";
    private String fittedTrackerHitCollectionName = "SVTFittedRawTrackerHits";
    private String trackerHitCollectionName = "TrackerHits";
    private String siClusterCollectionName = "StripClusterer_SiTrackerHitStrip1D";
    private String rotatedMCRelationsCollectionName = "RotatedHelicalTrackMCRelations";
    private String trackCollectionName = "MatchedTracks";
    private String trackerName = "Tracker";
    private Detector detector = null;
    IDDecoder dec;
    private IProfile1D peffFindable;
    private IProfile1D thetaeffFindable;
    private IProfile1D phieffFindable;
    private IProfile1D ctheffFindable;
    private IProfile1D d0effFindable;
    private IProfile1D z0effFindable;
    private IProfile1D peffElectrons;
    private IProfile1D thetaeffElectrons;
    private IProfile1D phieffElectrons;
    private IProfile1D ctheffElectrons;
    private IProfile1D d0effElectrons;
    private IProfile1D z0effElectrons;
    double beamP = 2.2;
    int nlayers = 12;
    int totelectrons = 0;
    double foundelectrons = 0;
    int findableelectrons = 0;
    int findableTracks = 0;
    double foundTracks = 0;
    private static final String nameStrip = "Tracker_TestRunModule_";
    private List<SiSensor> sensors;
    private boolean debugTrackEfficiency = false;

    public void setHelicalTrackHitCollectionName(String helicalTrackHitCollectionName) {
        this.helicalTrackHitCollectionName = helicalTrackHitCollectionName;
    }

    public void setTrackCollectionName(String trackCollectionName) {
        this.trackCollectionName = trackCollectionName;
    }

    public void setDebugTrackEfficiency(boolean debug) {
        this.debugTrackEfficiency = debug;
    }

    @Override
    protected void detectorChanged(Detector detector) {
        this.detector = detector;
        aida.tree().cd("/");
        IHistogramFactory hf = aida.histogramFactory();

        peffFindable = hf.createProfile1D("Findable Efficiency vs p", "", 20, 0., beamP);
        thetaeffFindable = hf.createProfile1D("Findable Efficiency vs theta", "", 20, 80, 100);
        phieffFindable = hf.createProfile1D("Findable Efficiency vs phi", "", 25, -0.25, 0.25);
        ctheffFindable = hf.createProfile1D("Findable Efficiency vs cos(theta)", "", 25, -0.25, 0.25);
        d0effFindable = hf.createProfile1D("Findable Efficiency vs d0", "", 50, -2., 2.);
        z0effFindable = hf.createProfile1D("Findable Efficiency vs z0", "", 50, -2., 2.);

        peffElectrons = hf.createProfile1D("Electrons Efficiency vs p", "", 20, 0., beamP);
        thetaeffElectrons = hf.createProfile1D("Electrons Efficiency vs theta", "", 20, 80, 100);
        phieffElectrons = hf.createProfile1D("Electrons Efficiency vs phi", "", 25, -0.25, 0.25);
        ctheffElectrons = hf.createProfile1D("Electrons Efficiency vs cos(theta)", "", 25, -0.25, 0.25);
        d0effElectrons = hf.createProfile1D("Electrons Efficiency vs d0", "", 20, -1., 1.);
        z0effElectrons = hf.createProfile1D("Electrons Efficiency vs z0", "", 20, -1., 1.);

        // Make a list of SiSensors in the SVT.
        sensors = this.detector.getSubdetector(trackerName).getDetectorElement().findDescendants(SiSensor.class);

        // Setup the occupancy plots.
        aida.tree().cd("/");
        for (int kk = 1; kk < 13; kk++) {
            IProfile1D clEffic = createLayerPlot("clusterEfficiency", kk, 50, 0, 25.);
        }
    }

    @Override
    public void process(EventHeader event) {

        aida.tree().cd("/");

        //make sure the required collections exist
        if (!event.hasCollection(RawTrackerHit.class, rawTrackerHitCollectionName))
            return;
        if (!event.hasCollection(FittedRawTrackerHit.class, fittedTrackerHitCollectionName))
            return;
        if (!event.hasCollection(Track.class, trackCollectionName))
            return;
        if (!event.hasCollection(LCRelation.class, rotatedMCRelationsCollectionName))
            return;
        if (!event.hasCollection(SiTrackerHitStrip1D.class, siClusterCollectionName))
            return;

        if (!event.hasCollection(SimTrackerHit.class, trackerHitCollectionName))
            return;
        //
        //get the b-field
        Hep3Vector IP = new BasicHep3Vector(0., 0., 1.);
        double bfield = event.getDetector().getFieldMap().getField(IP).y();
        //make some maps and relation tables        
        Map<Track, TrackAnalysis> tkanalMap = new HashMap<Track, TrackAnalysis>();
        RelationalTable hittomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        List<LCRelation> mcrelations = event.get(LCRelation.class, rotatedMCRelationsCollectionName);
        for (LCRelation relation : mcrelations) {
            if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                hittomc.add(relation.getFrom(), relation.getTo());
        }
        RelationalTable mcHittomcP = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        //  Get the collections of SimTrackerHits
        List<List<SimTrackerHit>> simcols = event.get(SimTrackerHit.class);
        //  Loop over the SimTrackerHits and fill in the relational table
        for (List<SimTrackerHit> simlist : simcols) {
            for (SimTrackerHit simhit : simlist) {
                if (simhit.getMCParticle() != null)
                    mcHittomcP.add(simhit, simhit.getMCParticle());
            }
        }
        RelationalTable trktomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, "SVTTrueHitRelations")) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, "SVTTrueHitRelations");
            for (LCRelation relation : trueHitRelations) {
                if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                    rawtomc.add(relation.getFrom(), relation.getTo());
            }
        }
        // make relational table for strip clusters to mc particle
        List<SiTrackerHitStrip1D> siClusters = event.get(SiTrackerHitStrip1D.class, siClusterCollectionName);
        RelationalTable clustertosimhit = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        for (SiTrackerHit cluster : siClusters) {
            List<RawTrackerHit> rawHits = cluster.getRawHits();
            for (RawTrackerHit rth : rawHits) {
                Set<SimTrackerHit> simTrackerHits = rawtomc.allFrom(rth);
                if (simTrackerHits != null)
                    for (SimTrackerHit simhit : simTrackerHits) {
                        clustertosimhit.add(cluster, simhit);
                    }
            }
        }
        //relational tables from mc particle to raw and fitted tracker hits
        RelationalTable fittomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        List<FittedRawTrackerHit> fittedTrackerHits = event.get(FittedRawTrackerHit.class, fittedTrackerHitCollectionName);
        for (FittedRawTrackerHit hit : fittedTrackerHits) {
            RawTrackerHit rth = hit.getRawTrackerHit();
            Set<SimTrackerHit> simTrackerHits = rawtomc.allFrom(rth);
            if (simTrackerHits != null)
                for (SimTrackerHit simhit : simTrackerHits) {
                    if (simhit.getMCParticle() != null)
                        fittomc.add(hit, simhit.getMCParticle());
                }
        }

        //  Instantiate the class that determines if a track is "findable"
        FindableTrack findable = new FindableTrack(event);

        List<Track> tracks = event.get(Track.class, trackCollectionName);
        aida.histogram1D("Tracks per Event").fill(tracks.size());
        for (Track trk : tracks) {
            TrackAnalysis tkanal = new TrackAnalysis(trk, hittomc);
            tkanalMap.put(trk, tkanal);
            MCParticle mcp = tkanal.getMCParticleNew();
            if (mcp != null)
                //  Create a map between the tracks found and the assigned MC particle
                trktomc.add(trk, tkanal.getMCParticleNew());
        }

        //  Now loop over all MC Particles
        List<MCParticle> mclist = event.getMCParticles();
        int _nchMCP = 0;
        int _nchMCPBar = 0;
        for (MCParticle mcp : mclist) {

            //  Calculate the pT and polar angle of the MC particle
            double px = mcp.getPX();
            double py = mcp.getPY();
            double pz = mcp.getPZ();
            double pt = Math.sqrt(px * px + py * py);
            double p = Math.sqrt(pt * pt + pz * pz);
            double cth = pz / p;
            double theta = 180. * Math.acos(cth) / Math.PI;
            double eta = -Math.log(Math.tan(Math.atan2(pt, pz) / 2));
            double phi = Math.atan2(py, px);
            //  Find the number of layers hit by this mc particle
//            System.out.println("MC pt=" + pt);
            int nhits = findable.LayersHit(mcp);
            boolean isFindable = findable.InnerTrackerIsFindable(mcp, nlayers - 2);

            //  Calculate the helix parameters for this MC particle
            HelixParamCalculator helix = new HelixParamCalculator(mcp, bfield);
            double d0 = helix.getDCA();
            double z0 = helix.getZ0();

            //  Check cases where we have multiple tracks associated with this MC particle
            Set<Track> trklist = trktomc.allTo(mcp);
            int ntrk = trklist.size();

//            Set<Track> trklistAxial = trktomcAxial.allTo(mcp);
//            int ntrkAxial = trklistAxial.size();
            if (mcp.getPDGID() == 622) {
                boolean bothreco = true;
                boolean bothfindable = true;
                //it's the A'...let's see if we found both tracks.
                List<MCParticle> daughters = mcp.getDaughters();
                for (MCParticle d : daughters) {
                    if (trktomc.allTo(d).isEmpty())
                        bothreco = false;
                    if (!findable.InnerTrackerIsFindable(d, nlayers - 2))
                        bothfindable = false;
                }
                double vtxWgt = 0;
                if (bothreco)
                    vtxWgt = 1.0;
//                VxEff.fill(mcp.getOriginX(), vtxWgt);
//                VyEff.fill(mcp.getOriginY(), vtxWgt);
//                VzEff.fill(mcp.getOriginZ(), vtxWgt);
                if (bothfindable) {
//                    VxEffFindable.fill(mcp.getOriginX(), vtxWgt);
//                    VyEffFindable.fill(mcp.getOriginY(), vtxWgt);
//                    VzEffFindable.fill(mcp.getOriginZ(), vtxWgt);
                }
            }

//            if (nhits == nlayers[0]) {
            if (isFindable) {
                _nchMCP++;
                findableTracks++;
                double wgt = 0.;
                if (ntrk > 0)
                    wgt = 1.;
                foundTracks += wgt;
                peffFindable.fill(p, wgt);
                phieffFindable.fill(phi, wgt);
                thetaeffFindable.fill(theta, wgt);
                ctheffFindable.fill(cth, wgt);
                d0effFindable.fill(d0, wgt);
                z0effFindable.fill(z0, wgt);

                if (wgt == 0) {
                    Set<SimTrackerHit> mchitlist = mcHittomcP.allTo(mcp);
                    Set<HelicalTrackCross> hitlist = hittomc.allTo(mcp);
                    Set<FittedRawTrackerHit> fitlist = fittomc.allTo(mcp);
                    if (debugTrackEfficiency) {
                        System.out.println("TrackMCEfficiencyMonitoring::  Missed a findable track with MC p = " + p);
                        if (!hasHTHInEachLayer(hitlist, fitlist))
                            System.out.println("This track failed becasue it's missing a helical track hit");
                    }
                }

            }
            if (mcp.getParents().size() == 1 && mcp.getParents().get(0).getPDGID() == 622) {
                totelectrons++;
//                    findableelectrons++;
                double wgt = 0.;
                if (ntrk > 0)
                    wgt = 1.;
                foundelectrons += wgt;
                peffElectrons.fill(p, wgt);
                phieffElectrons.fill(phi, wgt);
                thetaeffElectrons.fill(theta, wgt);
                ctheffElectrons.fill(cth, wgt);
                d0effElectrons.fill(d0, wgt);
                z0effElectrons.fill(z0, wgt);

                //               }
            }
        }
    }

    @Override
    public void fillEndOfRunPlots() {
    }

    @Override
    public void dumpDQMData() {
    }

    private IProfile1D getLayerPlot(String prefix, int layer) {
        return aida.profile1D(prefix + "_layer" + layer);
    }

    private IProfile1D createLayerPlot(String prefix, int layer, int nchan, double min, double max) {
        IProfile1D hist = aida.profile1D(prefix + "_layer" + layer, nchan, min, max);
        return hist;
    }

    private boolean hasHTHInEachLayer(Set<HelicalTrackCross> list, Set<FittedRawTrackerHit> fitlist) {
        for (int layer = 1; layer < nlayers - 2; layer += 2) {
            boolean hasThisLayer = false;
            for (HelicalTrackCross hit : list) {
                if (hit.Layer() == layer)
                    hasThisLayer = true;
            }
            if (!hasThisLayer) {
                System.out.println("Missing reconstructed hit in layer = " + layer);
                boolean hasFitHitSL1 = false;
                boolean hasFitHitSL2 = false;
                FittedRawTrackerHit fitSL1 = null;
                FittedRawTrackerHit fitSL2 = null;
                System.out.println("fitted hit list size = " + fitlist.size());
                for (FittedRawTrackerHit fit : fitlist) {
                    System.out.println("fitted hit layer number = " + fit.getRawTrackerHit().getLayerNumber());
                    if (fit.getRawTrackerHit().getLayerNumber() == layer) {
                        hasFitHitSL1 = true;
                        fitSL1 = fit;
                        System.out.println("Found a hit in SL1 with t0 = " + fitSL1.getT0() + "; amp = " + fitSL1.getAmp() + "; chi^2 = " + fitSL1.getShapeFitParameters().getChiSq() + "; strip = " + fitSL1.getRawTrackerHit().getCellID());
                    }
                    if (fit.getRawTrackerHit().getLayerNumber() == layer + 1) {
                        hasFitHitSL2 = true;
                        fitSL2 = fit;
                        System.out.println("Found a hit in SL2 with t0 = " + fitSL2.getT0() + "; amp = " + fitSL2.getAmp() + "; chi^2 = " + fitSL2.getShapeFitParameters().getChiSq() + "; strip = " + fitSL2.getRawTrackerHit().getCellID());

                    }
                }
                if (!hasFitHitSL1)
                    System.out.println("MISSING a hit in SL1!!!");
                if (!hasFitHitSL2)
                    System.out.println("MISSING a hit in SL2!!!");

                return false;
            }
        }
        return true;
    }

}
