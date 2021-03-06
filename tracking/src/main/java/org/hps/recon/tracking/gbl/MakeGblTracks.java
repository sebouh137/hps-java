package org.hps.recon.tracking.gbl;

import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.util.Pair;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.MultipleScattering;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.constants.Constants;
import org.lcsim.detector.ITransform3D;
import org.lcsim.detector.tracker.silicon.ChargeCarrier;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.detector.tracker.silicon.SiSensorElectrodes;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.base.BaseTrack;
import org.lcsim.event.base.BaseTrackState;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.fit.helicaltrack.HelicalTrackStrip;
import org.lcsim.fit.helicaltrack.HelixUtils;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.recon.tracking.digitization.sisim.TrackerHitType;

/**
 * Utilities that create track objects from fitted GBL trajectories.
 *
 * @author Per Hansson Adrian <phansson@slac.stanford.edu>
 *
 */
public class MakeGblTracks {

    private final static Logger LOGGER = Logger.getLogger(MakeGblTracks.class.getPackage().getName());
    static {
        LOGGER.setLevel(Level.WARNING);
    }
    
    
    private MakeGblTracks() {
    }

    public static void setDebug(boolean debug) {
        if (debug) {
            LOGGER.setLevel(Level.INFO);
        } else {
            LOGGER.setLevel(Level.OFF);
        }
    }

    /**
     * Create a new {@link BaseTrack} from a {@link FittedGblTrajectory}. 
     * @param fittedGblTrajectory
     * @param helicalTrackFit
     * @param hitsOnTrack
     * @param trackType
     * @param bfield
     * @return the new {@link BaseTrack} and the kinks along the {@link GblTrajectory} as a {@link Pair}.
     */
    public static Pair<Track, GBLKinkData> makeCorrectedTrack(FittedGblTrajectory fittedGblTrajectory, HelicalTrackFit helicalTrackFit, List<TrackerHit> hitsOnTrack, int trackType, double bfield) {
        //  Initialize the reference point to the origin
        double[] ref = new double[]{0., 0., 0.};

        //  Create a new SeedTrack
        BaseTrack trk = new BaseTrack();

        //  Add the hits to the track
        for (TrackerHit hit : hitsOnTrack) {
            trk.addHit(hit);
        }

        // Set base track parameters
        Pair<double[], SymmetricMatrix> correctedHelixParams = fittedGblTrajectory.getCorrectedPerigeeParameters(helicalTrackFit, FittedGblTrajectory.GBLPOINT.IP, bfield);
        trk.setTrackParameters(correctedHelixParams.getFirst(), bfield);// hack to set the track charge
        trk.getTrackStates().clear();

        // Set state at IP
        TrackState stateIP = new BaseTrackState(correctedHelixParams.getFirst(), ref, correctedHelixParams.getSecond().asPackedArray(true), TrackState.AtIP, bfield);
        trk.getTrackStates().add(stateIP);

        // Set state at last point on trajectory
        Pair<double[], SymmetricMatrix> correctedHelixParamsLast = fittedGblTrajectory.getCorrectedPerigeeParameters(helicalTrackFit, FittedGblTrajectory.GBLPOINT.LAST, bfield);
        TrackState stateLast = new BaseTrackState(correctedHelixParamsLast.getFirst(), ref, correctedHelixParamsLast.getSecond().asPackedArray(true), TrackState.AtLastHit, bfield);
        trk.getTrackStates().add(stateLast);

        // Extract kinks from trajectory
        GBLKinkData kinkData = fittedGblTrajectory.getKinks();

        // Set other info needed
        trk.setChisq(fittedGblTrajectory.get_chi2());
        trk.setNDF(fittedGblTrajectory.get_ndf());
        trk.setFitSuccess(true);
        trk.setRefPointIsDCA(true);
        trk.setTrackType(TrackType.setGBL(trackType, true));

        //  Add the track to the list of tracks
//            tracks.add(trk);
        LOGGER.fine(String.format("helix chi2 %f ndf %d gbl chi2 %f ndf %d\n", helicalTrackFit.chisqtot(), helicalTrackFit.ndf()[0] + helicalTrackFit.ndf()[1], trk.getChi2(), trk.getNDF()));
        return new Pair<Track, GBLKinkData>(trk, kinkData);
    }

    

    /**
     * Do a GBL fit to an arbitrary set of strip hits, with a starting value of
     * the helix parameters.
     *
     * @param helix Initial helix parameters. Only track parameters are used
     * (not covariance)
     * @param stripHits Strip hits to be used for the GBL fit. Does not need to
     * be in sorted order.
     * @param hth Stereo hits for the track's hit list (these are not used in
     * the GBL fit). Does not need to be in sorted order.
     * @param nIterations Number of times to iterate the GBL fit.
     * @param scattering Multiple scattering manager.
     * @param bfield B-field
     * @return The refitted track.
     */
    public static Pair<Track, GBLKinkData> refitTrack(HelicalTrackFit helix, Collection<TrackerHit> stripHits, Collection<TrackerHit> hth, int nIterations, int trackType, MultipleScattering scattering, double bfield) {
        List<TrackerHit> allHthList = TrackUtils.sortHits(hth);
        List<TrackerHit> sortedStripHits = TrackUtils.sortHits(stripHits);
        FittedGblTrajectory fit = doGBLFit(helix, sortedStripHits, scattering, bfield, 0);
        for (int i = 0; i < nIterations; i++) {
            Pair<Track, GBLKinkData> newTrack = makeCorrectedTrack(fit, helix, allHthList, trackType, bfield);
            helix = TrackUtils.getHTF(newTrack.getFirst());
            fit = doGBLFit(helix, sortedStripHits, scattering, bfield, 0);
        }
        Pair<Track, GBLKinkData> mergedTrack = makeCorrectedTrack(fit, helix, allHthList, trackType, bfield);
        return mergedTrack;
    }

    /**
     * Do a GBL fit to a list of {@link TrackerHit}. 
     * @param htf - seed fit
     * @param stripHits - list of {@link TrackerHit}.
     * @param _scattering - estimation of the multiple scattering {@link MultipleScattering}.
     * @param bfield - magnitude of B-field.
     * @param debug - debug flag.
     * @return the fitted GBL trajectory
     */
    public static FittedGblTrajectory doGBLFit(HelicalTrackFit htf, List<TrackerHit> stripHits, MultipleScattering _scattering, double bfield, int debug) {
        List<GBLStripClusterData> stripData = makeStripData(htf, stripHits, _scattering, bfield, debug);
        double bfac = Constants.fieldConversion * bfield;

        FittedGblTrajectory fit = HpsGblRefitter.fit(stripData, bfac, debug > 0);
        return fit;
    }

    /**
     * Create a list of {@link GBLStripClusterData} objects that can be used as input to the GBL fitter.
     * @param htf
     * @param stripHits
     * @param _scattering
     * @param _B
     * @param _debug
     * @return the list of GBL strip cluster data
     */
    public static List<GBLStripClusterData> makeStripData(HelicalTrackFit htf, List<TrackerHit> stripHits, MultipleScattering _scattering, double _B, int _debug) {
        List<GBLStripClusterData> stripClusterDataList = new ArrayList<GBLStripClusterData>();

        // Find scatter points along the path
        MultipleScattering.ScatterPoints scatters = _scattering.FindHPSScatterPoints(htf);

        if (_debug > 0) {
            System.out.printf("perPar covariance matrix\n%s\n", htf.covariance().toString());
        }

        for (TrackerHit stripHit : stripHits) {
            HelicalTrackStripGbl strip;
            if (stripHit instanceof SiTrackerHitStrip1D) {
                strip = new HelicalTrackStripGbl(makeDigiStrip((SiTrackerHitStrip1D) stripHit), true);
            } else {
                SiTrackerHitStrip1D newHit = new SiTrackerHitStrip1D(stripHit);
                strip = new HelicalTrackStripGbl(makeDigiStrip(newHit), true);
            }

            // find Millepede layer definition from DetectorElement
            HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) stripHit.getRawHits().get(0)).getDetectorElement();

            int millepedeId = sensor.getMillepedeId();

            if (_debug > 0) {
                System.out.printf("layer %d millepede %d (DE=\"%s\", origin %s) \n", strip.layer(), millepedeId, sensor.getName(), strip.origin().toString());
            }

            //Center of the sensor
            Hep3Vector origin = strip.origin();

            //Find intercept point with sensor in tracking frame
            Hep3Vector trkpos = TrackUtils.getHelixPlaneIntercept(htf, strip, Math.abs(_B));
            if (trkpos == null) {
                if (_debug > 0) {
                    System.out.println("Can't find track intercept; use sensor origin");
                }
                trkpos = strip.origin();
            }
            if (_debug > 0) {
                System.out.printf("trkpos at intercept [%.10f %.10f %.10f]\n", trkpos.x(), trkpos.y(), trkpos.z());
            }

            //GBLDATA
            GBLStripClusterData stripData = new GBLStripClusterData(millepedeId);
            //Add to output list
            stripClusterDataList.add(stripData);

            //path length to intercept
            double s = HelixUtils.PathToXPlane(htf, trkpos.x(), 0, 0).get(0);
            double s3D = s / Math.cos(Math.atan(htf.slope()));

            //GBLDATA
            stripData.setPath(s);
            stripData.setPath3D(s3D);

            //GBLDATA
            stripData.setU(strip.u());
            stripData.setV(strip.v());
            stripData.setW(strip.w());

            //Print track direction at intercept
            Hep3Vector tDir = HelixUtils.Direction(htf, s);
            double phi = htf.phi0() - s / htf.R();
            double lambda = Math.atan(htf.slope());

            //GBLDATA
            stripData.setTrackDir(tDir);
            stripData.setTrackPhi(phi);
            stripData.setTrackLambda(lambda);

            //Print residual in measurement system
            // start by find the distance vector between the center and the track position
            Hep3Vector vdiffTrk = VecOp.sub(trkpos, origin);

            // then find the rotation from tracking to measurement frame
            Hep3Matrix trkToStripRot = getTrackToStripRotation(sensor);

            // then rotate that vector into the measurement frame to get the predicted measurement position
            Hep3Vector trkpos_meas = VecOp.mult(trkToStripRot, vdiffTrk);

            //GBLDATA
            stripData.setMeas(strip.umeas());
            stripData.setTrackPos(trkpos_meas);
            stripData.setMeasErr(strip.du());

            if (_debug > 1) {
                System.out.printf("rotation matrix to meas frame\n%s\n", VecOp.toString(trkToStripRot));
                System.out.printf("tPosGlobal %s origin %s\n", trkpos.toString(), origin.toString());
                System.out.printf("tDiff %s\n", vdiffTrk.toString());
                System.out.printf("tPosMeas %s\n", trkpos_meas.toString());
            }

            if (_debug > 0) {
                System.out.printf("layer %d millePedeId %d uRes %.10f\n", strip.layer(), millepedeId, stripData.getMeas() - stripData.getTrackPos().x());
            }

            // find scattering angle
            MultipleScattering.ScatterPoint scatter = scatters.getScatterPoint(((RawTrackerHit) strip.getStrip().rawhits().get(0)).getDetectorElement());
            double scatAngle;

            if (scatter != null) {
                scatAngle = scatter.getScatterAngle().Angle();
            } else {
                if (_debug > 0) {
                    System.out.printf("WARNING cannot find scatter for detector %s with strip cluster at %s\n", ((RawTrackerHit) strip.getStrip().rawhits().get(0)).getDetectorElement().getName(), strip.origin().toString());
                }
                scatAngle = GblUtils.estimateScatter(sensor, htf, _scattering, _B);
            }

            //GBLDATA
            stripData.setScatterAngle(scatAngle);
        }
        return stripClusterDataList;
    }

    private static Hep3Matrix getTrackToStripRotation(SiSensor sensor) {
        // This function transforms the hit to the sensor coordinates

        // Transform from JLab frame to sensor frame (done through the RawTrackerHit)
        SiSensorElectrodes electrodes = sensor.getReadoutElectrodes(ChargeCarrier.HOLE);
        ITransform3D detToStrip = electrodes.getGlobalToLocal();
        // Get rotation matrix
        Hep3Matrix detToStripMatrix = detToStrip.getRotation().getRotationMatrix();
        // Transformation between the JLAB and tracking coordinate systems
        Hep3Matrix detToTrackMatrix = CoordinateTransformations.getMatrix();

        return VecOp.mult(detToStripMatrix, VecOp.inverse(detToTrackMatrix));
    }

    private static HelicalTrackStrip makeDigiStrip(SiTrackerHitStrip1D h) {
        SiTrackerHitStrip1D local = h.getTransformedHit(TrackerHitType.CoordinateSystem.SENSOR);
        SiTrackerHitStrip1D global = h.getTransformedHit(TrackerHitType.CoordinateSystem.GLOBAL);

        ITransform3D trans = local.getLocalToGlobal();
        Hep3Vector org = trans.transformed(new BasicHep3Vector(0., 0., 0.));
        Hep3Vector u = global.getMeasuredCoordinate();
        Hep3Vector v = global.getUnmeasuredCoordinate();

        //rotate to tracking frame
        Hep3Vector neworigin = CoordinateTransformations.transformVectorToTracking(org);
        Hep3Vector newu = CoordinateTransformations.transformVectorToTracking(u);
        Hep3Vector newv = CoordinateTransformations.transformVectorToTracking(v);

        double umeas = local.getPosition()[0];
        double vmin = VecOp.dot(local.getUnmeasuredCoordinate(), local.getHitSegment().getStartPoint());
        double vmax = VecOp.dot(local.getUnmeasuredCoordinate(), local.getHitSegment().getEndPoint());
        double du = Math.sqrt(local.getCovarianceAsMatrix().diagonal(0));

        //don't fill fields we don't use
//        IDetectorElement de = h.getSensor();
//        String det = getName(de);
//        int lyr = getLayer(de);
//        BarrelEndcapFlag be = getBarrelEndcapFlag(de);
        double dEdx = h.getdEdx();
        double time = h.getTime();
        List<RawTrackerHit> rawhits = h.getRawHits();
        HelicalTrackStrip strip = new HelicalTrackStrip(neworigin, newu, newv, umeas, du, vmin, vmax, dEdx, time, rawhits, null, -1, null);

        return strip;
    }
}
