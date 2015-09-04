package org.hps.recon.tracking.gbl;

import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.util.Pair;
import org.hps.recon.tracking.StrategyType;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.tracking.gbl.matrix.Matrix;
import org.hps.recon.tracking.gbl.matrix.SymMatrix;
import org.hps.recon.tracking.gbl.matrix.Vector;
import org.hps.util.BasicLogFormatter;
import org.lcsim.constants.Constants;
import org.lcsim.event.EventHeader;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.fit.helicaltrack.HelicalTrackHit;
import org.lcsim.fit.helicaltrack.HelixUtils;
import org.lcsim.lcio.LCIOConstants;
import org.lcsim.recon.tracking.seedtracker.SeedCandidate;
import org.lcsim.recon.tracking.seedtracker.SeedTrack;
import org.lcsim.util.log.LogUtil;

/**
 * A class that creates track objects from fitted GBL trajectories and adds them
 * into the event.
 *
 * @author Per Hansson Adrian <phansson@slac.stanford.edu>
 *
 */
public class MakeGblTracks {

    private String _TrkCollectionName = "GBLTracks";
    private static Logger logger = LogUtil.create(MakeGblTracks.class, new BasicLogFormatter());

    /**
     * Creates a new instance of MakeTracks.
     */
    public MakeGblTracks() {
        //logger = Logger.getLogger(getClass().getName());
        //logger.setUseParentHandlers(false);
        //Handler handler = new StreamHandler(System.out, new SimpleFormatter());
        //logger.addHandler(handler);
        logger.setLevel(Level.INFO);
//        try {
//            logger.addHandler(new FileHandler(MakeGblTracks.class.getSimpleName()+".log"));
//        } catch (SecurityException | IOException e) {
//            e.printStackTrace();
//        }

    }

    public void setDebug(boolean debug) {
        if (debug) {
            logger.setLevel(Level.INFO);
        } else {
            logger.setLevel(Level.OFF);
        }
    }

    /**
     * Process a Gbl track and store it into the event
     *
     * @param event event header
     * @param track Gbl trajectory
     * @param seed SeedTrack
     * @param bfield magnetic field (used to turn curvature into momentum)
     */
    public void Process(EventHeader event, List<FittedGblTrajectory> gblTrajectories, double bfield) {

        List<Track> tracks = new ArrayList<Track>();

        logger.info("adding " + gblTrajectories.size() + " of fitted GBL tracks to the event");

        for (FittedGblTrajectory fittedTraj : gblTrajectories) {

            //  Initialize the reference point to the origin
            double[] ref = new double[]{0., 0., 0.};
            SeedTrack seedTrack = (SeedTrack) fittedTraj.get_seed();
            SeedCandidate trackseed = seedTrack.getSeedCandidate();

            //  Create a new SeedTrack (SeedTrack extends BaseTrack)
            SeedTrack trk = new SeedTrack();

            //  Add the hits to the track
            for (HelicalTrackHit hit : trackseed.getHits()) {
                trk.addHit((TrackerHit) hit);
            }

            //  Retrieve the helix and save the relevant bits of helix info
            HelicalTrackFit helix = trackseed.getHelix();
            Pair<double[], SymmetricMatrix> correctedHelixParams = getGblCorrectedHelixParameters(helix, fittedTraj, bfield);
            trk.setTrackParameters(correctedHelixParams.getFirst(), bfield); // Sets first TrackState.
            trk.setCovarianceMatrix(correctedHelixParams.getSecond()); // Modifies first TrackState.
            trk.setChisq(fittedTraj.get_chi2());
            trk.setNDF(fittedTraj.get_ndf());

            //  Flag that the fit was successful and set the reference point
            trk.setFitSuccess(true);
            trk.setReferencePoint(ref); // Modifies first TrackState.
            trk.setRefPointIsDCA(true);

            //      Set the strategy used to find this track
            trk.setStratetgy(trackseed.getSeedStrategy());

            //  Set the SeedCandidate this track is based on
            trk.setSeedCandidate(trackseed);

            // Check if a StrategyType is associated with this strategy.
            // If it is, set the track type with the GBL flag set to true.
            // Otherwise, just move on and stick with the default value.
            StrategyType strategyType = StrategyType.getType(seedTrack.getType());
            if (strategyType != null) {
                trk.setTrackType(TrackType.getType(strategyType, true));
            }

            // Check the track - hook for plugging in external constraint
            //if ((_trackCheck != null) && (! _trackCheck.checkTrack(trk))) continue;
            //  Add the track to the list of tracks
            tracks.add((Track) trk);
            logger.info(String.format("helix chi2 %f ndf %d gbl chi2 %f ndf %d\n", helix.chisqtot(), helix.ndf()[0] + helix.ndf()[1], trk.getChi2(), trk.getNDF()));
            if (logger.getLevel().intValue() <= Level.INFO.intValue()) {
                for (int i = 0; i < 5; ++i) {
                    logger.info(String.format("param %d: %.10f -> %.10f    helix-gbl= %f", i, helix.parameters()[i], trk.getTrackParameter(i), helix.parameters()[i] - trk.getTrackParameter(i)));
                }
            }

        }

        logger.info("adding " + Integer.toString(tracks.size()) + " Gbl tracks to event with " + event.get(Track.class, "MatchedTracks").size() + " matched tracks");

        // Put the tracks back into the event and exit
        int flag = 1 << LCIOConstants.TRBIT_HITS;
        event.put(_TrkCollectionName, tracks, Track.class, flag);
    }

    /**
     * Compute the updated helix parameters and covariance matrix.
     *
     * @param helix - original seed track
     * @param traj - fitted GBL trajectory
     * @return corrected parameters.
     */
    private Pair<double[], SymmetricMatrix> getGblCorrectedHelixParameters(HelicalTrackFit helix, FittedGblTrajectory traj, double bfield) {

        // get seed helix parameters
//        double p = helix.p(Math.abs(bfield));
//        double q = Math.signum(helix.R());
//        double qOverP = q / p;
//        ClParams clParams = new ClParams(helix, bfield);
//        PerigeeParams perParams = new PerigeeParams(helix, bfield);
        double qOverP = helix.curvature() / (Constants.fieldConversion * Math.abs(bfield) * Math.sqrt(1 + Math.pow(helix.slope(), 2)));
        double d0 = -1.0 * helix.dca(); // correct for different sign convention of d0 in perigee frame
        double z0 = helix.z0();
        double phi0 = helix.phi0();
        double lambda = Math.atan(helix.slope());

//        System.out.println("clParams: " + clParams.getParams());
//        System.out.println("perParams: " + perParams.getParams());
//        System.out.format("converted params: qOverP %f, d0 %f, z0 %f, phi0 %f, lambda %f\n", qOverP, d0, z0, phi0, lambda);
        logger.info(String.format("original helix: d0=%f, z0=%f, omega=%f, tanlambda=%f, phi0=%f, p=%f", helix.dca(), helix.z0(), helix.curvature(), helix.slope(), helix.phi0(), helix.p(Math.abs(bfield))));
        logger.info("original helix covariance:\n" + helix.covariance());

        // get corrections from GBL fit
        Vector locPar = new Vector(5);
        SymMatrix locCov = new SymMatrix(5);
        traj.get_traj().getResults(1, locPar, locCov); // vertex point
//        locCov.print(10, 8);
        double qOverPCorr = locPar.get(FittedGblTrajectory.GBLPARIDX.QOVERP.getValue());
        double xTPrimeCorr = locPar.get(FittedGblTrajectory.GBLPARIDX.XTPRIME.getValue());
        double yTPrimeCorr = locPar.get(FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue());
        double xTCorr = locPar.get(FittedGblTrajectory.GBLPARIDX.XT.getValue());
        double yTCorr = locPar.get(FittedGblTrajectory.GBLPARIDX.YT.getValue());

        logger.info((helix.slope() > 0 ? "top: " : "bot ") + "qOverPCorr " + qOverPCorr + " xtPrimeCorr " + xTPrimeCorr + " yTPrimeCorr " + yTPrimeCorr + " xTCorr " + xTCorr + " yTCorr " + yTCorr);

        // calculate new d0 and z0
        Hep3Matrix perToClPrj = traj.get_track_data().getPrjPerToCl();
        Hep3Matrix clToPerPrj = VecOp.inverse(perToClPrj);
        Hep3Vector corrPer = VecOp.mult(clToPerPrj, new BasicHep3Vector(xTCorr, yTCorr, 0.0));

        //d0
        double d0_corr = corrPer.y();
        double dca_gbl = -1.0 * (d0 + d0_corr);

        //z0
        double z0_corr = corrPer.z();
        double z0_gbl = z0 + z0_corr;

        //calculate new phi0
        double phi0_gbl = phi0 + xTPrimeCorr;

        //calculate new slope
        double lambda_gbl = lambda + yTPrimeCorr;
        double slope_gbl = Math.tan(lambda_gbl);

        // calculate new curvature
        double qOverP_gbl = qOverP + qOverPCorr;
//        double pt_gbl = (1.0 / qOverP_gbl) * Math.cos(lambda_gbl);
//        double C_gbl = Constants.fieldConversion * Math.abs(bfield) / pt_gbl;
        double C_gbl = Constants.fieldConversion * Math.abs(bfield) * qOverP_gbl / Math.cos(lambda_gbl);

        logger.info("qOverP=" + qOverP + " qOverPCorr=" + qOverPCorr + " qOverP_gbl=" + qOverP_gbl + " ==> pGbl=" + 1.0 / qOverP_gbl + " C_gbl=" + C_gbl);

        logger.info(String.format("corrected helix: d0=%f, z0=%f, omega=%f, tanlambda=%f, phi0=%f, p=%f", dca_gbl, z0_gbl, C_gbl, slope_gbl, phi0_gbl, Math.abs(1 / qOverP_gbl)));

        // Strandlie, Wittek, NIMA 566, 2006
        Matrix covariance_gbl = new Matrix(5, 5);
        //helpers
        double Bz = -Constants.fieldConversion * Math.abs(bfield); // TODO sign convention and should it be it scaled from Telsa?
        double p = Math.abs(1 / qOverP_gbl);
        double q = Math.signum(qOverP_gbl);
        double tanLambda = Math.tan(lambda_gbl);
        double cosLambda = Math.cos(lambda_gbl);
        Hep3Vector B = new BasicHep3Vector(0, 0, 1); // TODO sign convention?
        Hep3Vector H = VecOp.mult(1 / bfield, B);
        Hep3Vector T = HelixUtils.Direction(helix, 0.);
        Hep3Vector HcrossT = VecOp.cross(H, T);
        double alpha = HcrossT.magnitude(); // this should be Bvec cross TrackDir/|B|
        double Q = Math.abs(bfield) * q / p;
        Hep3Vector Z = new BasicHep3Vector(0, 0, 1);
        Hep3Vector J = VecOp.mult(1. / VecOp.cross(T, Z).magnitude(), VecOp.cross(T, Z));
        Hep3Vector K = Z;
        Hep3Vector U = VecOp.mult(-1, J);
        Hep3Vector V = VecOp.cross(T, U);
        Hep3Vector I = VecOp.cross(J, K);
        Hep3Vector N = VecOp.mult(1 / alpha, VecOp.cross(H, T));
        double UdotI = VecOp.dot(U, I);
        double NdotV = VecOp.dot(N, V);
        double NdotU = VecOp.dot(N, U);
        double TdotI = VecOp.dot(T, I);
        double VdotI = VecOp.dot(V, I);
        double VdotK = VecOp.dot(V, K);
        covariance_gbl.set(HelicalTrackFit.dcaIndex, FittedGblTrajectory.GBLPARIDX.XT.getValue(), VdotK / TdotI);
        covariance_gbl.set(HelicalTrackFit.phi0Index, FittedGblTrajectory.GBLPARIDX.XTPRIME.getValue(), 1);
        covariance_gbl.set(HelicalTrackFit.phi0Index, FittedGblTrajectory.GBLPARIDX.XT.getValue(), -alpha * Q * UdotI * NdotU / (cosLambda * TdotI));
        covariance_gbl.set(HelicalTrackFit.phi0Index, FittedGblTrajectory.GBLPARIDX.YT.getValue(), -alpha * Q * VdotI * NdotU / (cosLambda * TdotI));
        covariance_gbl.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.QOVERP.getValue(), -1 * Bz / cosLambda);
//        covariance_gbl.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.XTPRIME.getValue(), 0);
        covariance_gbl.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue(), -1 * q * Bz * tanLambda / (p * cosLambda));
        covariance_gbl.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.XT.getValue(), q * Bz * alpha * Q * tanLambda * UdotI * NdotV / (p * cosLambda * TdotI));
        covariance_gbl.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.YT.getValue(), q * Bz * alpha * Q * tanLambda * VdotI * NdotV / (p * cosLambda * TdotI));
        covariance_gbl.set(HelicalTrackFit.z0Index, FittedGblTrajectory.GBLPARIDX.YT.getValue(), -1 / TdotI);
        covariance_gbl.set(HelicalTrackFit.slopeIndex, FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue(), -1);
        covariance_gbl.set(HelicalTrackFit.slopeIndex, FittedGblTrajectory.GBLPARIDX.XT.getValue(), alpha * Q * UdotI * NdotV / TdotI);
        covariance_gbl.set(HelicalTrackFit.slopeIndex, FittedGblTrajectory.GBLPARIDX.YT.getValue(), alpha * Q * VdotI * NdotV / TdotI);

//        System.out.println(clToPerPrj);

//        covariance_gbl.print(15, 13);
////        System.out.println(-alpha * Q * UdotI * NdotU / (cosLambda * TdotI));
//        System.out.println(-alpha * Q * VdotI * NdotU / (cosLambda * TdotI));
//        System.out.format("%f %f %f %f %f %f\n", -alpha, Q, VdotI, NdotU, cosLambda, TdotI);
//        System.out.format("%f %f %f %f %f %f\n", -Math.cos(lambda)/Math.abs(bfield), Q, perToClPrj.e(1, 0), perToClPrj.e(0, 1), Math.cos(lambda_gbl), perToClPrj.e(2, 0));
//
////        System.out.println(q * Bz * alpha * Q * tanLambda * UdotI * NdotV / (p * cosLambda * TdotI));
//        System.out.println(q * Bz * alpha * Q * tanLambda * VdotI * NdotV / (p * cosLambda * TdotI));
////        System.out.println(alpha * Q * UdotI * NdotV / TdotI);
//        System.out.println(alpha * Q * VdotI * NdotV / TdotI);
        // Sho's magic
        Matrix jacobian = new Matrix(5, 5);
        jacobian.set(HelicalTrackFit.dcaIndex, FittedGblTrajectory.GBLPARIDX.XT.getValue(), -clToPerPrj.e(1, 0));
        jacobian.set(HelicalTrackFit.dcaIndex, FittedGblTrajectory.GBLPARIDX.YT.getValue(), -clToPerPrj.e(1, 1));
        jacobian.set(HelicalTrackFit.phi0Index, FittedGblTrajectory.GBLPARIDX.XTPRIME.getValue(), 1.0);
        jacobian.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.QOVERP.getValue(), Constants.fieldConversion * Math.abs(bfield) / Math.cos(lambda_gbl));
        jacobian.set(HelicalTrackFit.curvatureIndex, FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue(), Constants.fieldConversion * Math.abs(bfield) * qOverP_gbl * Math.tan(lambda_gbl) / Math.cos(lambda_gbl));
        jacobian.set(HelicalTrackFit.z0Index, FittedGblTrajectory.GBLPARIDX.XT.getValue(), clToPerPrj.e(2, 0));
        jacobian.set(HelicalTrackFit.z0Index, FittedGblTrajectory.GBLPARIDX.YT.getValue(), clToPerPrj.e(2, 1));
        jacobian.set(HelicalTrackFit.slopeIndex, FittedGblTrajectory.GBLPARIDX.YTPRIME.getValue(), Math.pow(Math.cos(lambda_gbl), -2.0));

//        jacobian.print(15, 13);
//        System.out.println(-clToPerPrj.e(1, 1));
//        System.out.println(clToPerPrj.e(2, 0));
        Matrix helixCovariance = jacobian.times(locCov.times(jacobian.transpose()));
        SymmetricMatrix cov = new SymmetricMatrix(5);
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                if (i >= j) {
                    cov.setElement(i, j, helixCovariance.get(i, j));
                }
            }
        }
        logger.info("corrected helix covariance:\n" + cov);

        double parameters_gbl[] = new double[5];
        parameters_gbl[HelicalTrackFit.dcaIndex] = dca_gbl;
        parameters_gbl[HelicalTrackFit.phi0Index] = phi0_gbl;
        parameters_gbl[HelicalTrackFit.curvatureIndex] = C_gbl;
        parameters_gbl[HelicalTrackFit.z0Index] = z0_gbl;
        parameters_gbl[HelicalTrackFit.slopeIndex] = slope_gbl;

        return new Pair<double[], SymmetricMatrix>(parameters_gbl, cov);
    }

    public void setTrkCollectionName(String name) {
        _TrkCollectionName = name;
    }

}
