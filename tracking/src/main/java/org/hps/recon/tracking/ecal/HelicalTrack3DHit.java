package org.hps.recon.tracking.ecal;

import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.BasicHep3Vector;
import java.util.List;
import org.lcsim.fit.helicaltrack.HelicalTrackHit;
import org.lcsim.geometry.subdetector.BarrelEndcapFlag;
import org.lcsim.fit.helicaltrack.HitUtils; 

/**
 * Encapsulate 3D hit info needed by HelicalTrackFitter to handle ECal hits (or otherwise non-layer hits).
 * This class is explicitly for HPS where the length of the
 * sensors are (mostly) along the detector
 * y-dimension ( == HelicalTrackFit x-dimension);
 * @author Bradley Yale <btu29@wildcats.unh.edu>
 * Copied/Modified from org.lcsim.recon.tracking.helicaltrack.HelicalTrack3DHit.java
 */

public class HelicalTrack3DHit extends HelicalTrackHit {
    private double _dz;
    private static int _type = 1;
    
    public HelicalTrack3DHit(Hep3Vector pos, SymmetricMatrix cov, double dEdx, double time,
            List rawhits, String detname, int layer, BarrelEndcapFlag beflag) {
        super(pos, cov, dEdx, time, _type, rawhits, detname, layer, beflag);
        _dz = Math.sqrt(cov.e(2, 2));
        if (! (_dz>_eps)) //_eps inherited from HelicalTrackHit
            _dz = _eps; 
    }
    
    /**
     * Create a HelicalTrack3DHit from a TrackerHit.
     * @param hit TrackerHit associated with this hit
     * @param beflag BarrelEndcapFlag for this hit
     */
//    public HelicalTrack3DHit(TrackerHit hit, BarrelEndcapFlag beflag) {
//        super(hit);
//        super.setBarrelEndcapFlag(beflag);
//        _dz = Math.sqrt(hit.getCovMatrix()[5]);
//    }
    
    /**
     * Create a HelicalTrack3DHit from a TrackerHit overriding the hit position
     * and covariance matrix in the TrackerHit.
     * @param hit TrackerHit associated with this hit
     * @param pos hit position
     * @param cov covariance matrix
     * @param beflag BarrelEndcapFlag for this hit
     */
//    public HelicalTrack3DHit(TrackerHit hit, Hep3Vector pos, SymmetricMatrix cov, BarrelEndcapFlag beflag) {
//        super(hit, pos, cov);
//        super.setBarrelEndcapFlag(beflag);
//        _dz = Math.sqrt(cov.e(2,2));
//    }
    
    /**
     * Create a HelicalTrack3DHit from scratch without reference to an existing
     * TrackerHit.  The BarrelEndcapFlag defaults to barrel, but can be set to
     * an endcap disk using the parent method setBarrelEndcapFlag.
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     * @param drphi uncertainty in the r*phi coordinate
     * @param dz uncertainty in the z coordinate
     */
    @Deprecated
    public HelicalTrack3DHit(double x, double y, double z, double drphi, double dz){
        this(new BasicHep3Vector(x, y, z), HitUtils.PixelCov(x, y, drphi, dz), 0., 0., null,
                "Unknown", 0 , BarrelEndcapFlag.BARREL);
        _dz = dz;
    }
    
    /**
     * Return the uncertainty in the z coordinate.
     * @return uncertainty in the z coordinate
     */
    public double dz() {
        if (super.BarrelEndcapFlag() == BarrelEndcapFlag.BARREL) return _dz;
        else throw new RuntimeException("z coordinate uncertainty undefined for a disk hit");
    }
}
    
