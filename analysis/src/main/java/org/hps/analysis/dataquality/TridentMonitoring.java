package org.hps.analysis.dataquality;

import hep.aida.IAnalysisFactory;
import hep.aida.IFitFactory;
import hep.aida.IFitResult;
import hep.aida.IFitter;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.particle.ReconParticleDriver;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.Vertex;
import org.lcsim.geometry.Detector;

/**
 * DQM driver V0 particles (i.e. e+e- pars) plots things like number of vertex
 * position an mass
 *
 * @author mgraham on May 14, 2014
 *
 */
public class TridentMonitoring extends DataQualityMonitor {

    private double ebeam = 1.05;
    private final BasicHep3Matrix beamAxisRotation = new BasicHep3Matrix();
    private static final int nCuts = 9;
    private final String[] cutNames = {"Trk Quality",
        "V0 Quality",
        "V0 Vertex",
        "Timing",
        "Tracking",
        "Cluster",
        "Event",
        "Front Hits",
        "Isolation"};
    private int firstVertexingCut = 0;
    private static final int PASS = 0;
    private static final int FAIL = 1;

    private final String finalStateParticlesColName = "FinalStateParticles";
    private final String unconstrainedV0CandidatesColName = "UnconstrainedV0Candidates";
//    private final String beamConV0CandidatesColName = "BeamspotConstrainedV0Candidates";
//    private final String targetV0ConCandidatesColName = "TargetConstrainedV0Candidates";
//    private final String trackListName = "MatchedTracks";
    private String[] fpQuantNames = {"nV0_per_Event", "avg_BSCon_mass", "avg_BSCon_Vx", "avg_BSCon_Vy", "avg_BSCon_Vz", "sig_BSCon_Vx", "sig_BSCon_Vy", "sig_BSCon_Vz", "avg_BSCon_Chi2"};

    private final String plotDir = "TridentMonitoring/";

    private IHistogram2D triTrackTime2D;
    private IHistogram1D triTrackTimeDiff;
    private IHistogram2D triMassMomentum;
    private IHistogram2D triZVsMomentum;
    private IHistogram2D triTrackMomentum2D;
    private IHistogram2D triPyEleVsPyPos;
    private IHistogram2D triPxEleVsPxPos;
    private IHistogram1D triDeltaP;
    private IHistogram1D triSumP;
    private IHistogram1D triMass;
    private IHistogram2D triZVsMass;
    private IHistogram1D triX;
    private IHistogram1D triY;
    private IHistogram1D triZ;
    private IHistogram2D triZY;
    private IHistogram2D triXY;
    private IHistogram1D triPx;
    private IHistogram1D triPy;
    private IHistogram1D triPz;
    private IHistogram2D triPxPy;
    private IHistogram1D triU;
    private IHistogram1D triV;

    private IHistogram2D triRadTrackTime2D;
    private IHistogram1D triRadTrackTimeDiff;
    private IHistogram2D triRadMassMomentum;
    private IHistogram2D triRadZVsMomentum;
    private IHistogram2D triRadTrackMomentum2D;
    private IHistogram2D triRadPyEleVsPyPos;
    private IHistogram2D triRadPxEleVsPxPos;
    private IHistogram1D triRadDeltaP;
    private IHistogram1D triRadSumP;
    private IHistogram1D triRadMass;
    private IHistogram2D triRadZVsMass;
    private IHistogram1D triRadX;
    private IHistogram1D triRadY;
    private IHistogram1D triRadZ;
    private IHistogram2D triRadZY;
    private IHistogram2D triRadXY;
    private IHistogram1D triRadPx;
    private IHistogram1D triRadPy;
    private IHistogram1D triRadPz;
    private IHistogram2D triRadPxPy;
    private IHistogram1D triRadU;
    private IHistogram1D triRadV;

    private IHistogram2D vertTrackTime2D;
    private IHistogram1D vertTrackTimeDiff;
    private IHistogram2D vertMassMomentum;
    private IHistogram2D vertZVsMomentum;
    private IHistogram2D vertTrackMomentum2D;
    private IHistogram2D vertPyEleVsPyPos;
    private IHistogram2D vertPxEleVsPxPos;
    private IHistogram1D vertDeltaP;
    private IHistogram1D vertSumP;
    private IHistogram1D vertMass;
    private IHistogram2D vertZVsMass;
    private IHistogram1D vertX;
    private IHistogram1D vertY;
    private IHistogram1D vertZ;
    private IHistogram2D vertZY;
    private IHistogram2D vertXY;
    private IHistogram1D vertPx;
    private IHistogram1D vertPy;
    private IHistogram1D vertPz;
    private IHistogram2D vertPxPy;
    private IHistogram1D vertU;
    private IHistogram1D vertV;

    private IHistogram2D vertRadTrackTime2D;
    private IHistogram1D vertRadTrackTimeDiff;
    private IHistogram2D vertRadMassMomentum;
    private IHistogram2D vertRadZVsMomentum;
    private IHistogram2D vertRadTrackMomentum2D;
    private IHistogram2D vertRadPyEleVsPyPos;
    private IHistogram2D vertRadPxEleVsPxPos;
    private IHistogram1D vertRadDeltaP;
    private IHistogram1D vertRadSumP;
    private IHistogram1D vertRadMass;
    private IHistogram2D vertRadZVsMass;
    private IHistogram1D vertRadX;
    private IHistogram1D vertRadY;
    private IHistogram1D vertRadZ;
    private IHistogram2D vertRadZY;
    private IHistogram2D vertRadXY;
    private IHistogram1D vertRadPx;
    private IHistogram1D vertRadPy;
    private IHistogram1D vertRadPz;
    private IHistogram2D vertRadPxPy;
    private IHistogram1D vertRadU;
    private IHistogram1D vertRadV;

    private IHistogram1D nTriCand;
    private IHistogram1D nVtxCand;
//    IHistogram1D vertexW;
//    IHistogram2D vertexVZ;

    private IHistogram1D maxTrkChi2;
    private IHistogram2D zVsMaxTrkChi2;
    private IHistogram1D v0Chi2;
    private IHistogram2D zVsV0Chi2;
    private IHistogram1D trackTimeDiff;
    private IHistogram2D zVsTrackTimeDiff;
    private IHistogram1D hitTimeStdDev;
    private IHistogram2D zVsHitTimeStdDev;
    private IHistogram1D eventTrkCount;
    private IHistogram1D eventPosCount;
    private IHistogram2D zVsEventTrkCount;
    private IHistogram2D zVsEventPosCount;
    private IHistogram1D l1Iso;
    private IHistogram2D zVsL1Iso;

    private final IHistogram1D[][] cutVertexMass = new IHistogram1D[nCuts][2];
    private final IHistogram1D[][] cutVertexZ = new IHistogram1D[nCuts][2];
    private final IHistogram2D[][] cutVertexZVsMass = new IHistogram2D[nCuts][2];

    //clean up event first
    private int nTrkMax = 5;
    private int nPosMax = 1;

    private double maxChi2SeedTrack = 7.0;
    private double maxChi2GBLTrack = 15.0;
    private double maxVertChi2 = 7.0;

    //v0 cuts   
    private double v0PzMax = 1.25 * ebeam;//GeV 
    private double v0PzMin = 0.1;// GeV
    private double v0PyMax = 0.04;//GeV absolute value
    private double v0PxMax = 0.04;//GeV absolute value
    private double v0VzMax = 50.0;// mm from target...someday make mass dependent
    private double v0VyMax = 1.0;// mm from target...someday make mass dependent
    private double v0VxMax = 2.0;// mm from target...someday make mass dependent
    //  track quality cuts
    private double beamPCut = 0.85;
    private double minPCut = 0.05;
//    private double trkPyMax = 0.2;
//    private double trkPxMax = 0.2;
    private double radCut = 0.8 * ebeam;
    private double trkTimeDiff = 5.0;
    private double clusterTimeDiffCut = 2.5;

    private double l1IsoMin = 1.0;
//cluster matching
//    private boolean reqCluster = false;
//    private int nClustMax = 3;
//    private double eneLossFactor = 0.7; //average E/p roughly
//    private double eneOverPCut = 0.3; //|(E/p)_meas - (E/p)_mean|<eneOverPCut

//counters
    private float nEvents = 0;
    private float nRecoV0 = 0;
    private final float[] nPassCut = new float[nCuts];

    public void setEbeam(double ebeam) {
        this.ebeam = ebeam;
    }

    @Override
    protected void detectorChanged(Detector detector) {
        System.out.println("TridendMonitoring::detectorChanged  Setting up the plotter");
        beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);

        aida.tree().cd("/");
        String trkType = "SeedTrack/";
        if (isGBL) {
            trkType = "GBLTrack/";
        }
        /*  V0 Quantities   */
        /*  Mass, vertex, chi^2 of fit */
        /* beamspot constrained */
//        IHistogram1D nV0 = aida.histogram1D(plotDir +  triggerType + "/"+"Number of V0 per event", 10, 0, 10);
//        IHistogram1D bsconMass = aida.histogram1D(plotDir +  triggerType + "/"+"BS Constrained Mass (GeV)", 100, 0, 0.200);
//        IHistogram1D bsconVx = aida.histogram1D(plotDir +  triggerType + "/"+"BS Constrained Vx (mm)", 50, -1, 1);
//        IHistogram1D bsconVy = aida.histogram1D(plotDir +  triggerType + "/"+"BS Constrained Vy (mm)", 50, -1, 1);
//        IHistogram1D bsconVz = aida.histogram1D(plotDir +  triggerType + "/"+"BS Constrained Vz (mm)", 50, -10, 10);
//        IHistogram1D bsconChi2 = aida.histogram1D(plotDir +  triggerType + "/"+"BS Constrained Chi2", 25, 0, 25);
//        /* target constrained */
//        IHistogram1D tarconMass = aida.histogram1D(plotDir +  triggerType + "/"+"Target Constrained Mass (GeV)", 100, 0, 0.200);
//        IHistogram1D tarconVx = aida.histogram1D(plotDir +  triggerType + "/"+ triggerType + "/"+"Target Constrained Vx (mm)", 50, -1, 1);
//        IHistogram1D tarconVy = aida.histogram1D(plotDir +  triggerType + "/"+ triggerType + "/"+"Target Constrained Vy (mm)", 50, -1, 1);
//        IHistogram1D tarconVz = aida.histogram1D(plotDir +  triggerType + "/"+ triggerType + "/"+"Target Constrained Vz (mm)", 50, -10, 10);
//        IHistogram1D tarconChi2 = aida.histogram1D(plotDir +  triggerType + "/"+ triggerType + "/"+"Target Constrained Chi2", 25, 0, 25);
        nTriCand = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Number of Trident Candidates", 5, 0, 4);

        triTrackTimeDiff = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Track time difference", 100, -10, 10);
        triTrackTime2D = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Track time vs. track time", 100, -10, 10, 100, -10, 10);

        triTrackMomentum2D = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Positron vs. electron momentum", 100, 0, v0PzMax, 100, 0, v0PzMax);
        triDeltaP = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Positron - electron momentum", 100, -1., 1.0);
        triSumP = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Positron + electron momentum", 100, v0PzMin, v0PzMax);
        triPyEleVsPyPos = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Py(e) vs Py(p)", 50, -0.04, 0.04, 50, -0.04, 0.04);
        triPxEleVsPxPos = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Px(e) vs Px(p)", 50, -0.04, 0.04, 50, -0.04, 0.04);

        triMassMomentum = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Vertex mass vs. vertex momentum", 100, v0PzMin, v0PzMax, 100, 0, 0.1);
        triZVsMomentum = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Z vs. vertex momentum", 100, v0PzMin, v0PzMax, 100, -v0VzMax, v0VzMax);
        triMass = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex mass", 100, 0, 0.11);
        triZVsMass = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Z vs. mass", 100, 0, 0.11, 100, -v0VzMax, v0VzMax);
        triX = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex X", 100, -v0VxMax, v0VxMax);
        triY = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Y", 100, -v0VyMax, v0VyMax);
        triZ = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Z", 100, -v0VzMax, v0VzMax);
        triXY = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Y vs. X", 100, -v0VxMax, v0VxMax, 100, -v0VyMax, v0VyMax);
        triZY = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Z vs. Y", 100, -v0VyMax, v0VyMax, 100, -v0VzMax, v0VzMax);
        triPx = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Px", 100, -v0PxMax, v0PxMax);
        triPy = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Py", 100, -v0PyMax, v0PyMax);
        triPz = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Pz", 100, v0PzMin, v0PzMax);
        triPxPy = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Py vs. Px", 100, -v0PxMax, v0PxMax, 100, -v0PyMax, v0PyMax);
        triU = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Px over Ptot", 100, -0.1, 0.1);
        triV = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Trident: Vertex Py over Ptot", 100, -0.1, 0.1);

        triRadTrackTimeDiff = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Track time difference", 100, -10, 10);
        triRadTrackTime2D = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Track time vs. track time", 100, -10, 10, 100, -10, 10);

        triRadTrackMomentum2D = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Positron vs. electron momentum", 100, 0, v0PzMax, 100, 0, v0PzMax);
        triRadDeltaP = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Positron - electron momentum", 100, -1., 1.0);
        triRadSumP = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Positron + electron momentum", 100, v0PzMin, v0PzMax);
        triRadPyEleVsPyPos = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Py(e) vs Py(p)", 50, -0.04, 0.04, 50, -0.04, 0.04);
        triRadPxEleVsPxPos = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Px(e) vs Px(p)", 50, -0.04, 0.04, 50, -0.04, 0.04);

        triRadMassMomentum = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex mass vs. vertex momentum", 100, v0PzMin, v0PzMax, 100, 0, 0.1);
        triRadZVsMomentum = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Z vs. vertex momentum", 100, v0PzMin, v0PzMax, 100, -v0VzMax, v0VzMax);
        triRadMass = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex mass", 100, 0, 0.11);
        triRadZVsMass = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Z vs. mass", 100, 0, 0.11, 100, -v0VzMax, v0VzMax);
        triRadX = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex X", 100, -v0VxMax, v0VxMax);
        triRadY = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Y", 100, -v0VyMax, v0VyMax);
        triRadZ = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Z", 100, -v0VzMax, v0VzMax);
        triRadXY = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Y vs. X", 100, -v0VxMax, v0VxMax, 100, -v0VyMax, v0VyMax);
        triRadZY = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Z vs. Y", 100, -v0VyMax, v0VyMax, 100, -v0VzMax, v0VzMax);
        triRadPx = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Px", 100, -v0PxMax, v0PxMax);
        triRadPy = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Py", 100, -v0PyMax, v0PyMax);
        triRadPz = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Pz", 100, v0PzMin, v0PzMax);
        triRadPxPy = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Py vs. Px", 100, -v0PxMax, v0PxMax, 100, -v0PyMax, v0PyMax);
        triRadU = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Px over Ptot", 100, -0.1, 0.1);
        triRadV = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative trident: Vertex Py over Ptot", 100, -0.1, 0.1);

        vertTrackTimeDiff = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Track time difference", 100, -10, 10);
        vertTrackTime2D = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Track time vs. track time", 100, -10, 10, 100, -10, 10);

        vertTrackMomentum2D = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Positron vs. electron momentum", 100, 0, v0PzMax, 100, 0, v0PzMax);
        vertDeltaP = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Positron - electron momentum", 100, -1., 1.0);
        vertSumP = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Positron + electron momentum", 100, v0PzMin, v0PzMax);
        vertPyEleVsPyPos = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Py(e) vs Py(p)", 50, -0.04, 0.04, 50, -0.04, 0.04);
        vertPxEleVsPxPos = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Px(e) vs Px(p)", 50, -0.04, 0.04, 50, -0.04, 0.04);

        vertMassMomentum = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex mass vs. vertex momentum", 100, v0PzMin, v0PzMax, 100, 0, 0.1);
        vertZVsMomentum = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Z vs. vertex momentum", 100, v0PzMin, v0PzMax, 100, -v0VzMax, v0VzMax);
        vertMass = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex mass", 100, 0, 0.11);
        vertZVsMass = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Z vs. mass", 100, 0, 0.11, 100, -v0VzMax, v0VzMax);
        vertX = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex X", 100, -v0VxMax, v0VxMax);
        vertY = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Y", 100, -v0VyMax, v0VyMax);
        vertZ = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Z", 100, -v0VzMax, v0VzMax);
        vertXY = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Y vs. X", 100, -v0VxMax, v0VxMax, 100, -v0VyMax, v0VyMax);
        vertZY = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Z vs. Y", 100, -v0VyMax, v0VyMax, 100, -v0VzMax, v0VzMax);
        vertPx = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Px", 100, -v0PxMax, v0PxMax);
        vertPy = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Py", 100, -v0PyMax, v0PyMax);
        vertPz = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Pz", 100, v0PzMin, v0PzMax);
        vertPxPy = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Py vs. Px", 100, -v0PxMax, v0PxMax, 100, -v0PyMax, v0PyMax);
        vertU = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Px over Ptot", 100, -0.1, 0.1);
        vertV = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Vertex: Vertex Py over Ptot", 100, -0.1, 0.1);

        vertRadTrackTimeDiff = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Track time difference", 100, -10, 10);
        vertRadTrackTime2D = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Track time vs. track time", 100, -10, 10, 100, -10, 10);

        vertRadTrackMomentum2D = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Positron vs. electron momentum", 100, 0, v0PzMax, 100, 0, v0PzMax);
        vertRadDeltaP = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Positron - electron momentum", 100, -1., 1.0);
        vertRadSumP = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Positron + electron momentum", 100, v0PzMin, v0PzMax);
        vertRadPyEleVsPyPos = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Py(e) vs Py(p)", 50, -0.04, 0.04, 50, -0.04, 0.04);
        vertRadPxEleVsPxPos = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Px(e) vs Px(p)", 50, -0.04, 0.04, 50, -0.04, 0.04);

        vertRadMassMomentum = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex mass vs. vertex momentum", 100, v0PzMin, v0PzMax, 100, 0, 0.1);
        vertRadZVsMomentum = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Z vs. vertex momentum", 100, v0PzMin, v0PzMax, 100, -v0VzMax, v0VzMax);
        vertRadMass = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex mass", 100, 0, 0.11);
        vertRadZVsMass = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Z vs. mass", 100, 0, 0.11, 100, -v0VzMax, v0VzMax);
        vertRadX = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex X", 100, -v0VxMax, v0VxMax);
        vertRadY = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Y", 100, -v0VyMax, v0VyMax);
        vertRadZ = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Z", 100, -v0VzMax, v0VzMax);
        vertRadXY = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Y vs. X", 100, -v0VxMax, v0VxMax, 100, -v0VyMax, v0VyMax);
        vertRadZY = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Z vs. Y", 100, -v0VyMax, v0VyMax, 100, -v0VzMax, v0VzMax);
        vertRadPx = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Px", 100, -v0PxMax, v0PxMax);
        vertRadPy = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Py", 100, -v0PyMax, v0PyMax);
        vertRadPz = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Pz", 100, v0PzMin, v0PzMax);
        vertRadPxPy = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Py vs. Px", 100, -v0PxMax, v0PxMax, 100, -v0PyMax, v0PyMax);
        vertRadU = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Px over Ptot", 100, -0.1, 0.1);
        vertRadV = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Radiative vertex: Vertex Py over Ptot", 100, -0.1, 0.1);

        nVtxCand = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Number of Vertexing Candidates", 5, 0, 4);

        maxTrkChi2 = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Cut: Trk Chi2", 50, 0.0, 50.0);
        zVsMaxTrkChi2 = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Cut: Vz vs Trk Chi2", 50, 0.0, 50.0, 50, -v0VzMax, v0VzMax);

        v0Chi2 = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Cut: V0 Chi2", 50, 0.0, 25.0);
        zVsV0Chi2 = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Cut: Vz vs V0 Chi2", 50, 0.0, 25.0, 50, -v0VzMax, v0VzMax);

        trackTimeDiff = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Cut: Trk Time Diff", 50, 0.0, 10.0);
        hitTimeStdDev = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Cut: Hit Time Std Dev", 50, 0.0, 10.0);
        zVsTrackTimeDiff = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Cut: Vz vs Trk Time Diff", 50, 0.0, 10.0, 50, -v0VzMax, v0VzMax);
        zVsHitTimeStdDev = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Cut: Vz vs Hit Time Std Dev", 50, 0.0, 10.0, 50, -v0VzMax, v0VzMax);

        eventTrkCount = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Cut: Num Tracks", 10, 0.5, 10.5);
        eventPosCount = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Cut: Num Positrons", 5, 0.5, 5.5);
        zVsEventTrkCount = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Cut: Vz vs Num Tracks", 10, 0.5, 10.5, 50, -v0VzMax, v0VzMax);
        zVsEventPosCount = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Cut: Vz vs Num Positrons", 5, 0.5, 5.5, 50, -v0VzMax, v0VzMax);

        l1Iso = aida.histogram1D(plotDir + trkType + triggerType + "/" + "Cut: L1 Isolation", 50, 0.0, 5.0);
        zVsL1Iso = aida.histogram2D(plotDir + trkType + triggerType + "/" + "Cut: Vz vs L1 Isolation", 50, 0.0, 5.0, 50, -v0VzMax, v0VzMax);

        for (int i = 0; i < nCuts; i++) {
            for (int pass = 0; pass < 2; pass++) {
                cutVertexZ[i][pass] = aida.histogram1D(String.format("%s%s%s/cut %d: %s/%s: Vertex Z position (mm)", plotDir, trkType, triggerType, i, cutNames[i], pass == PASS ? "pass" : "fail"),
                        100, -v0VzMax, v0VzMax);
                cutVertexMass[i][pass] = aida.histogram1D(String.format("%s%s%s/cut %d: %s/%s: Vertex mass (GeV)", plotDir, trkType, triggerType, i, cutNames[i], pass == PASS ? "pass" : "fail"),
                        100, 0, 0.1 * ebeam);
                cutVertexZVsMass[i][pass] = aida.histogram2D(String.format("%s%s%s/cut %d: %s/%s: Vertex Z vs. mass", plotDir, trkType, triggerType, i, cutNames[i], pass == PASS ? "pass" : "fail"),
                        100, 0, 0.1 * ebeam, 100, -v0VzMax, v0VzMax);
            }
        }
    }

    @Override
    public void process(EventHeader event) {
        /*  make sure everything is there */
        if (!event.hasCollection(ReconstructedParticle.class, finalStateParticlesColName)) {
            return;
        }
        if (!event.hasCollection(ReconstructedParticle.class, unconstrainedV0CandidatesColName)) {
            return;
        }
//        if (!event.hasCollection(ReconstructedParticle.class, beamConV0CandidatesColName)) {
//            return;
//        }
//        if (!event.hasCollection(ReconstructedParticle.class, targetV0ConCandidatesColName)) {
//            return;
//        }

        //check to see if this event is from the correct trigger (or "all");
        if (!matchTrigger(event)) {
            return;
        }

        nEvents++;

        int nV0 = 0;
        List<ReconstructedParticle> unConstrainedV0List = event.get(ReconstructedParticle.class, unconstrainedV0CandidatesColName);
        for (ReconstructedParticle v0 : unConstrainedV0List) {
            if (isGBL == TrackType.isGBL(v0.getType())) {
                nV0++;
            }
        }
        nRecoV0 += nV0;

        RelationalTable hitToStrips = TrackUtils.getHitToStripsTable(event);
        RelationalTable hitToRotated = TrackUtils.getHitToRotatedTable(event);

        List<ReconstructedParticle> fspList = event.get(ReconstructedParticle.class, finalStateParticlesColName);
        int npos = 0;
        int ntrk = 0;
        for (ReconstructedParticle fsp : fspList) {
            if (isGBL != TrackType.isGBL(fsp.getType())) {
                continue;
            }
            if (fsp.getCharge() != 0) {
                ntrk++;
            }
            if (fsp.getCharge() > 0) {
                npos++;
            }
        }

        List<ReconstructedParticle> candidateList = new ArrayList<>();
        List<ReconstructedParticle> vertCandidateList = new ArrayList<>();
        for (ReconstructedParticle uncV0 : unConstrainedV0List) {
            if (isGBL != TrackType.isGBL(uncV0.getType())) {
                continue;
            }
            Vertex uncVert = uncV0.getStartVertex();
//  v0 & vertex-quality cuts
            Hep3Vector v0MomRot = VecOp.mult(beamAxisRotation, uncV0.getMomentum());
            Hep3Vector v0Vtx = VecOp.mult(beamAxisRotation, uncVert.getPosition());

            List<Track> tracks = new ArrayList<Track>();
            ReconstructedParticle electron = uncV0.getParticles().get(ReconParticleDriver.ELECTRON);
            ReconstructedParticle positron = uncV0.getParticles().get(ReconParticleDriver.POSITRON);
            if (electron.getCharge() != -1 || positron.getCharge() != 1) {
                throw new RuntimeException("incorrect charge on v0 daughters");
            }
            tracks.add(electron.getTracks().get(0));
            tracks.add(positron.getTracks().get(0));
            if (tracks.size() != 2) {
                throw new RuntimeException("expected two tracks in vertex, got " + tracks.size());
            }
            List<Double> trackTimes = new ArrayList<Double>();
            List<Double> hitTimes = new ArrayList<Double>();
            double mean = 0;
            for (Track track : tracks) {
                trackTimes.add(TrackUtils.getTrackTime(track, hitToStrips, hitToRotated));
                for (TrackerHit hit : TrackUtils.getStripHits(track, hitToStrips, hitToRotated)) {
                    mean += hit.getTime();
                    hitTimes.add(hit.getTime());
                }
            }
            mean /= hitTimes.size();
            double stdDev = 0;
            for (Double time : hitTimes) {
                stdDev += Math.pow(time - mean, 2);
            }
            stdDev /= (hitTimes.size() - 1);
            stdDev = Math.sqrt(stdDev);

            Double[] eleIso = TrackUtils.getIsolations(electron.getTracks().get(0), hitToStrips, hitToRotated);
            Double[] posIso = TrackUtils.getIsolations(positron.getTracks().get(0), hitToStrips, hitToRotated);
            double minL1Iso = -9999;
            if (eleIso[0] != null && posIso[0] != null) {
                double eleL1Iso = Math.min(Math.abs(eleIso[0]), Math.abs(eleIso[1]));
                double posL1Iso = Math.min(Math.abs(posIso[0]), Math.abs(posIso[1]));
                minL1Iso = Math.min(eleL1Iso, posL1Iso);
            }

            //start applying cuts
            int cutNum = 0;
            boolean trackQualityCut = Math.max(tracks.get(0).getChi2(), tracks.get(1).getChi2()) < (isGBL ? maxChi2GBLTrack : maxChi2SeedTrack);
            maxTrkChi2.fill(Math.max(tracks.get(0).getChi2(), tracks.get(1).getChi2()));
            zVsMaxTrkChi2.fill(Math.max(tracks.get(0).getChi2(), tracks.get(1).getChi2()), v0Vtx.z());
            if (!trackQualityCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            boolean v0QualityCut = uncVert.getChi2() < maxVertChi2;
            v0Chi2.fill(uncVert.getChi2());
            zVsV0Chi2.fill(uncVert.getChi2(), v0Vtx.z());
            if (!v0QualityCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            boolean vertexMomentumCut = v0MomRot.z() < v0PzMax && v0MomRot.z() > v0PzMin && Math.abs(v0MomRot.x()) < v0PxMax && Math.abs(v0MomRot.y()) < v0PyMax;
            boolean vertexPositionCut = Math.abs(v0Vtx.x()) < v0VxMax && Math.abs(v0Vtx.y()) < v0VyMax && Math.abs(v0Vtx.z()) < v0VzMax;
            if (!vertexMomentumCut || !vertexPositionCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            boolean trackTimeDiffCut = Math.abs(trackTimes.get(0) - trackTimes.get(1)) < trkTimeDiff;
            trackTimeDiff.fill(Math.abs(trackTimes.get(0) - trackTimes.get(1)));
            hitTimeStdDev.fill(stdDev);
            zVsTrackTimeDiff.fill(Math.abs(trackTimes.get(0) - trackTimes.get(1)), v0Vtx.z());
            zVsHitTimeStdDev.fill(stdDev, v0Vtx.z());
            if (!trackTimeDiffCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            boolean topBottomCut = electron.getMomentum().y() * positron.getMomentum().y() < 0;
            boolean pMinCut = electron.getMomentum().magnitude() > minPCut && positron.getMomentum().magnitude() > minPCut;
            boolean pMaxCut = electron.getMomentum().magnitude() < beamPCut && positron.getMomentum().magnitude() < beamPCut;
            if (!topBottomCut || !pMaxCut || !pMinCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            boolean clusterMatchCut = !electron.getClusters().isEmpty() && !positron.getClusters().isEmpty();
            boolean clusterTimeCut = clusterMatchCut && Math.abs(ClusterUtilities.getSeedHitTime(electron.getClusters().get(0)) - ClusterUtilities.getSeedHitTime(positron.getClusters().get(0))) < clusterTimeDiffCut;
            if (!clusterMatchCut || !clusterTimeCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            boolean eventTrkCountCut = ntrk >= 2 && ntrk <= nTrkMax;
            boolean eventPosCountCut = npos >= 1 && npos <= nPosMax;
            eventTrkCount.fill(ntrk);
            eventPosCount.fill(npos);
            zVsEventTrkCount.fill(ntrk, v0Vtx.z());
            zVsEventPosCount.fill(npos, v0Vtx.z());
            if (!eventTrkCountCut || !eventPosCountCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            candidateList.add(uncV0);
            firstVertexingCut = cutNum;

            boolean frontHitsCut = eleIso[0] != null && posIso[0] != null && eleIso[2] != null && posIso[2] != null;
            if (!frontHitsCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            l1Iso.fill(minL1Iso);
            zVsL1Iso.fill(minL1Iso, v0Vtx.z());
            boolean isoCut = minL1Iso > l1IsoMin;
            if (!isoCut) {
                cutVertexZ[cutNum][FAIL].fill(v0Vtx.z());
                cutVertexMass[cutNum][FAIL].fill(uncV0.getMass());
                cutVertexZVsMass[cutNum][FAIL].fill(uncV0.getMass(), v0Vtx.z());
                continue;
            }
            cutVertexZ[cutNum][PASS].fill(v0Vtx.z());
            cutVertexMass[cutNum][PASS].fill(uncV0.getMass());
            cutVertexZVsMass[cutNum][PASS].fill(uncV0.getMass(), v0Vtx.z());
            nPassCut[cutNum]++;
            cutNum++;

            vertCandidateList.add(uncV0);
        }

        nTriCand.fill(candidateList.size());
        nVtxCand.fill(vertCandidateList.size());

        if (!candidateList.isEmpty()) {
            // pick the best candidate...for now just pick a random one. 
            ReconstructedParticle bestCandidate = candidateList.get((int) (Math.random() * candidateList.size()));

            //fill some stuff: 
            ReconstructedParticle electron = bestCandidate.getParticles().get(ReconParticleDriver.ELECTRON);
            ReconstructedParticle positron = bestCandidate.getParticles().get(ReconParticleDriver.POSITRON);
            if (electron.getCharge() != -1 || positron.getCharge() != 1) {
                throw new RuntimeException("vertex needs e+ and e- but is missing one or both");
            }

            double tEle = TrackUtils.getTrackTime(electron.getTracks().get(0), hitToStrips, hitToRotated);
            double tPos = TrackUtils.getTrackTime(positron.getTracks().get(0), hitToStrips, hitToRotated);
            Hep3Vector pBestV0Rot = VecOp.mult(beamAxisRotation, bestCandidate.getMomentum());
            Hep3Vector pEleRot = VecOp.mult(beamAxisRotation, electron.getMomentum());
            Hep3Vector pPosRot = VecOp.mult(beamAxisRotation, positron.getMomentum());
            Hep3Vector v0Vtx = VecOp.mult(beamAxisRotation, bestCandidate.getStartVertex().getPosition());

            triTrackTime2D.fill(tEle, tPos);
            triTrackTimeDiff.fill(tEle - tPos);
            triZVsMomentum.fill(bestCandidate.getMomentum().magnitude(), v0Vtx.z());
            triMassMomentum.fill(bestCandidate.getMomentum().magnitude(), bestCandidate.getMass());
            triTrackMomentum2D.fill(electron.getMomentum().magnitude(), positron.getMomentum().magnitude());
            triPyEleVsPyPos.fill(pEleRot.y(), pPosRot.y());
            triPxEleVsPxPos.fill(pEleRot.x(), pPosRot.x());
            triSumP.fill(bestCandidate.getMomentum().magnitude());
            triDeltaP.fill(positron.getMomentum().magnitude() - electron.getMomentum().magnitude());

            triPxPy.fill(pBestV0Rot.x(), pBestV0Rot.y());
            triMass.fill(bestCandidate.getMass());
            triZVsMass.fill(bestCandidate.getMass(), v0Vtx.z());
            triX.fill(v0Vtx.x());
            triY.fill(v0Vtx.y());
            triZ.fill(v0Vtx.z());
            triPx.fill(pBestV0Rot.x());
            triPy.fill(pBestV0Rot.y());
            triPz.fill(pBestV0Rot.z());
            triU.fill(pBestV0Rot.x() / pBestV0Rot.magnitude());
            triV.fill(pBestV0Rot.y() / pBestV0Rot.magnitude());
            triXY.fill(v0Vtx.x(), v0Vtx.y());
            triZY.fill(v0Vtx.y(), v0Vtx.z());
            if (bestCandidate.getMomentum().magnitude() > radCut) {
                triRadTrackTime2D.fill(tEle, tPos);
                triRadTrackTimeDiff.fill(tEle - tPos);
                triRadZVsMomentum.fill(bestCandidate.getMomentum().magnitude(), v0Vtx.z());
                triRadMassMomentum.fill(bestCandidate.getMomentum().magnitude(), bestCandidate.getMass());
                triRadTrackMomentum2D.fill(electron.getMomentum().magnitude(), positron.getMomentum().magnitude());
                triRadPyEleVsPyPos.fill(pEleRot.y(), pPosRot.y());
                triRadPxEleVsPxPos.fill(pEleRot.x(), pPosRot.x());
                triRadSumP.fill(bestCandidate.getMomentum().magnitude());
                triRadDeltaP.fill(positron.getMomentum().magnitude() - electron.getMomentum().magnitude());

                triRadPxPy.fill(pBestV0Rot.x(), pBestV0Rot.y());
                triRadMass.fill(bestCandidate.getMass());
                triRadZVsMass.fill(bestCandidate.getMass(), v0Vtx.z());
                triRadX.fill(v0Vtx.x());
                triRadY.fill(v0Vtx.y());
                triRadZ.fill(v0Vtx.z());
                triRadPx.fill(pBestV0Rot.x());
                triRadPy.fill(pBestV0Rot.y());
                triRadPz.fill(pBestV0Rot.z());
                triRadU.fill(pBestV0Rot.x() / pBestV0Rot.magnitude());
                triRadV.fill(pBestV0Rot.y() / pBestV0Rot.magnitude());
                triRadXY.fill(v0Vtx.x(), v0Vtx.y());
                triRadZY.fill(v0Vtx.y(), v0Vtx.z());
            }
        }

        if (!vertCandidateList.isEmpty()) {
            // pick the best candidate...for now just pick a random one. 
            ReconstructedParticle bestCandidate = vertCandidateList.get((int) (Math.random() * vertCandidateList.size()));

            //fill some stuff: 
            ReconstructedParticle electron = bestCandidate.getParticles().get(ReconParticleDriver.ELECTRON);
            ReconstructedParticle positron = bestCandidate.getParticles().get(ReconParticleDriver.POSITRON);
            if (electron.getCharge() != -1 || positron.getCharge() != 1) {
                throw new RuntimeException("vertex needs e+ and e- but is missing one or both");
            }

            double tEle = TrackUtils.getTrackTime(electron.getTracks().get(0), hitToStrips, hitToRotated);
            double tPos = TrackUtils.getTrackTime(positron.getTracks().get(0), hitToStrips, hitToRotated);
            Hep3Vector pBestV0Rot = VecOp.mult(beamAxisRotation, bestCandidate.getMomentum());
            Hep3Vector pEleRot = VecOp.mult(beamAxisRotation, electron.getMomentum());
            Hep3Vector pPosRot = VecOp.mult(beamAxisRotation, positron.getMomentum());
            Hep3Vector v0Vtx = VecOp.mult(beamAxisRotation, bestCandidate.getStartVertex().getPosition());

            vertTrackTime2D.fill(tEle, tPos);
            vertTrackTimeDiff.fill(tEle - tPos);
            vertZVsMomentum.fill(bestCandidate.getMomentum().magnitude(), v0Vtx.z());
            vertMassMomentum.fill(bestCandidate.getMomentum().magnitude(), bestCandidate.getMass());
            vertTrackMomentum2D.fill(electron.getMomentum().magnitude(), positron.getMomentum().magnitude());
            vertPyEleVsPyPos.fill(pEleRot.y(), pPosRot.y());
            vertPxEleVsPxPos.fill(pEleRot.x(), pPosRot.x());
            vertSumP.fill(bestCandidate.getMomentum().magnitude());
            vertDeltaP.fill(positron.getMomentum().magnitude() - electron.getMomentum().magnitude());

            vertPxPy.fill(pBestV0Rot.x(), pBestV0Rot.y());
            vertMass.fill(bestCandidate.getMass());
            vertZVsMass.fill(bestCandidate.getMass(), v0Vtx.z());
            vertX.fill(v0Vtx.x());
            vertY.fill(v0Vtx.y());
            vertZ.fill(v0Vtx.z());
            vertPx.fill(pBestV0Rot.x());
            vertPy.fill(pBestV0Rot.y());
            vertPz.fill(pBestV0Rot.z());
            vertU.fill(pBestV0Rot.x() / pBestV0Rot.magnitude());
            vertV.fill(pBestV0Rot.y() / pBestV0Rot.magnitude());
            vertXY.fill(v0Vtx.x(), v0Vtx.y());
            vertZY.fill(v0Vtx.y(), v0Vtx.z());
            if (bestCandidate.getMomentum().magnitude() > radCut) {
                vertRadTrackTime2D.fill(tEle, tPos);
                vertRadTrackTimeDiff.fill(tEle - tPos);
                vertRadZVsMomentum.fill(bestCandidate.getMomentum().magnitude(), v0Vtx.z());
                vertRadMassMomentum.fill(bestCandidate.getMomentum().magnitude(), bestCandidate.getMass());
                vertRadTrackMomentum2D.fill(electron.getMomentum().magnitude(), positron.getMomentum().magnitude());
                vertRadPyEleVsPyPos.fill(pEleRot.y(), pPosRot.y());
                vertRadPxEleVsPxPos.fill(pEleRot.x(), pPosRot.x());
                vertRadSumP.fill(bestCandidate.getMomentum().magnitude());
                vertRadDeltaP.fill(positron.getMomentum().magnitude() - electron.getMomentum().magnitude());

                vertRadPxPy.fill(pBestV0Rot.x(), pBestV0Rot.y());
                vertRadMass.fill(bestCandidate.getMass());
                vertRadZVsMass.fill(bestCandidate.getMass(), v0Vtx.z());
                vertRadX.fill(v0Vtx.x());
                vertRadY.fill(v0Vtx.y());
                vertRadZ.fill(v0Vtx.z());
                vertRadPx.fill(pBestV0Rot.x());
                vertRadPy.fill(pBestV0Rot.y());
                vertRadPz.fill(pBestV0Rot.z());
                vertRadU.fill(pBestV0Rot.x() / pBestV0Rot.magnitude());
                vertRadV.fill(pBestV0Rot.y() / pBestV0Rot.magnitude());
                vertRadXY.fill(v0Vtx.x(), v0Vtx.y());
                vertRadZY.fill(v0Vtx.y(), v0Vtx.z());
            }
        }
    }

    @Override
    public void printDQMData() {
        System.out.println("TridendMonitoring::printDQMData");
        for (Entry<String, Double> entry : monitoredQuantityMap.entrySet()) {
            System.out.println(entry.getKey() + " = " + entry.getValue());
        }
        System.out.println("*******************************");

        System.out.println("TridendMonitoring::Tridend Selection Summary: " + (isGBL ? "GBLTrack" : "SeedTrack"));

        System.out.println("\t\t\tTridend Selection Summary");
        System.out.println("******************************************************************************************");
        System.out.format("Number of      V0:\t%8.0f\t%8.6f\t%8.6f\t%8.6f\n", nRecoV0, nRecoV0 / nRecoV0, nRecoV0 / nRecoV0, nRecoV0 / nEvents);

        for (int i = 0; i < nCuts; i++) {
            if (i == firstVertexingCut) {
                System.out.println("******************************************************************************************");
                System.out.println("\t\t\tVertex Selection Summary");
                System.out.println("******************************************************************************************");
            }
            System.out.format("%-12s Cuts:\t%8.0f\t%8.6f\t%8.6f\t%8.6f\n", cutNames[i], nPassCut[i], nPassCut[i] / (i == 0 ? nRecoV0 : nPassCut[i == 0 ? 0 : (i - 1)]), nPassCut[i] / nRecoV0, nPassCut[i] / nEvents);
        }
        System.out.println("******************************************************************************************");
    }

    /**
     * Calculate the averages here and fill the map
     */
    @Override
    public void calculateEndOfRunQuantities() {

        IAnalysisFactory analysisFactory = IAnalysisFactory.create();
        IFitFactory fitFactory = analysisFactory.createFitFactory();
        IFitter fitter = fitFactory.createFitter("chi2");

    }

    @Override
    public void printDQMStrings() {
        for (int i = 0; i < 9; i++)//TODO:  do this in a smarter way...loop over the map
        {
            System.out.println("ALTER TABLE dqm ADD " + fpQuantNames[i] + " double;");
        }
    }

    IFitResult fitVertexPosition(IHistogram1D h1d, IFitter fitter, double[] init, String range
    ) {
        return fitter.fit(h1d, "g+p1", init, range);
    }

}
