package org.hps.recon.tracking;

import java.util.List;

import org.lcsim.detector.converter.compact.subdetector.HpsTracker2;
import org.lcsim.detector.converter.compact.subdetector.SvtStereoLayer;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.EventHeader;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.base.BaseTrack;
import org.lcsim.event.base.BaseTrackState;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.util.Driver;

/**
 * Driver used to incrementally extrapolate the track to all SVT layers, 
 * using a fieldmap and the positions from the compact, and save the TrackState at each layer.
 * 
 * @author <a href="mailto:btu29@wildcats.unh.edu">Bradley Yale</a>
 * @author <a href="mailto:omoreno@slac.stanford.edu">Omar Moreno</a>
 */
public final class TrackStateDriver extends Driver {

    /** 
     * Name of the constants denoting the positions of the SVT layers and Ecal face 
     * in the compact description.
     */
	
 //   private static final String L1t_POSITION_CONSTANT_AXIAL_NAME = "module_L1t_halfmodule_axial_position";
 //   private static final String L1b_POSITION_CONSTANT_AXIAL_NAME = "module_L1b_halfmodule_axial_position"; 
 //   private static final String L2t_POSITION_CONSTANT_AXIAL_NAME = "module_L2t_halfmodule_axial_position";
 //   private static final String L2b_POSITION_CONSTANT_AXIAL_NAME = "module_L2b_halfmodule_axial_position";
 //   private static final String L3t_POSITION_CONSTANT_AXIAL_NAME = "module_L3t_halfmodule_axial_position";
 //   private static final String L3b_POSITION_CONSTANT_AXIAL_NAME = "module_L3b_halfmodule_axial_position";
 //   private static final String L4t_POSITION_CONSTANT_AXIAL_NAME = "module_L4t_halfmodule_axial_position";
 //   private static final String L4b_POSITION_CONSTANT_AXIAL_NAME = "module_L4b_halfmodule_axial_position";
 //   private static final String L5t_POSITION_CONSTANT_AXIAL_NAME = "module_L5t_halfmodule_axial_position";
 //   private static final String L5b_POSITION_CONSTANT_AXIAL_NAME = "module_L5b_halfmodule_axial_position";
 //   private static final String L6t_POSITION_CONSTANT_AXIAL_NAME = "module_L6t_halfmodule_axial_position";
 //   private static final String L6b_POSITION_CONSTANT_AXIAL_NAME = "module_L6b_halfmodule_axial_position";
    
 //   private static final String L1t_POSITION_CONSTANT_STEREO_NAME = "module_L1t_halfmodule_stereo_position";
 //   private static final String L1b_POSITION_CONSTANT_STEREO_NAME = "module_L1b_halfmodule_stereo_position"; 
 //   private static final String L2t_POSITION_CONSTANT_STEREO_NAME = "module_L2t_halfmodule_stereo_position";
 //   private static final String L2b_POSITION_CONSTANT_STEREO_NAME = "module_L2b_halfmodule_stereo_position";
 //   private static final String L3t_POSITION_CONSTANT_STEREO_NAME = "module_L3t_halfmodule_stereo_position";
 //   private static final String L3b_POSITION_CONSTANT_STEREO_NAME = "module_L3b_halfmodule_stereo_position";
 //   private static final String L4t_POSITION_CONSTANT_STEREO_NAME = "module_L4t_halfmodule_stereo_position";
 //   private static final String L4b_POSITION_CONSTANT_STEREO_NAME = "module_L4b_halfmodule_stereo_position";
 //   private static final String L5t_POSITION_CONSTANT_STEREO_NAME = "module_L5t_halfmodule_stereo_position";
 //   private static final String L5b_POSITION_CONSTANT_STEREO_NAME = "module_L5b_halfmodule_stereo_position";
 //   private static final String L6t_POSITION_CONSTANT_STEREO_NAME = "module_L6t_halfmodule_stereo_position";
 //   private static final String L6b_POSITION_CONSTANT_STEREO_NAME = "module_L6b_halfmodule_stereo_position";
    
    private static final String ECAL_POSITION_CONSTANT_NAME = "ecal_dface";

    /** Name of the SVT subdetector volume. */
    private static final String SUBDETECTOR_NAME = "Tracker";
    
    /** The B field map */
    FieldMap bFieldMap = null;
    
    /** The magnitude of the B field used.  */
    private double bField = 0.24; // Tesla

    private double ecalPosition = 0; // mm
   
    /** Z position to start extrapolation from */
    private double extStartPos = 700; // mm

    /** The extrapolation step size */ 
    private double stepSize = 5.0; // mm
    
    /** Top/Bottom SVT layer 1/2 z position */
    private double topLayer1Z = 0;
    private double botLayer1Z = 0;
    private double topLayer2Z = 0;
    private double botLayer2Z = 0;
    
    /** Top/Bottom SVT layer 3/4 z position */
    private double topLayer3Z = 0;
    private double botLayer3Z = 0;
    private double topLayer4Z = 0;
    private double botLayer4Z = 0;
    
    /** Top/Bottom SVT layer 5/6 z position */
    private double topLayer5Z = 0;
    private double botLayer5Z = 0;
    private double topLayer6Z = 0;
    private double botLayer6Z = 0;
    
    /** Top/Bottom SVT layer 7/8 z position */
    private double topLayer7Z = 0;
    private double botLayer7Z = 0;
    private double topLayer8Z = 0;
    private double botLayer8Z = 0;
    
    /** Top/Bottom SVT layer 9/10 z position */
    private double topLayer9Z = 0;
    private double botLayer9Z = 0;
    private double topLayer10Z = 0;
    private double botLayer10Z = 0;
    
    /** Top/Bottom SVT layer 11/12 z position */
    private double topLayer11Z = 0;
    private double botLayer11Z = 0;
    private double topLayer12Z = 0;
    private double botLayer12Z = 0;
    
    
    /** Name of the collection of tracks to apply corrections to. */
    private String gblTrackCollectionName = "GBLTracks";
   
    /** Name of the collection of seed tracks. */
    private String seedTrackCollectionName = "MatchedTracks";

    /** Default constructor */
    public TrackStateDriver() {}
      
    @Override
    protected void detectorChanged(Detector detector) {
      
        // Get the field map from the detector object
        bFieldMap = detector.getFieldMap(); 
        
        // Get the B-field from the geometry description 
        bField = TrackUtils.getBField(detector).magnitude();
       
        ecalPosition = detector.getConstants().get(ECAL_POSITION_CONSTANT_NAME).getValue();
    
        // Get the stereo layers from the geometry and build the stereo
        // layer maps
        List<SvtStereoLayer> stereoLayers 
            = ((HpsTracker2) detector.getSubdetector(SUBDETECTOR_NAME).getDetectorElement()).getStereoPairs();

        // Loop through all of the stereo layers and collect their sensor positions.
        // This will be used to set the track states at those layers.
        
        for (SvtStereoLayer stereoLayer : stereoLayers) { 
            if (stereoLayer.getLayerNumber() > 2) continue;
            
            // 1st stereo pair (layers 1 & 2)
            HpsSiSensor axialSensor1 = stereoLayer.getAxialSensor();
            HpsSiSensor stereoSensor1 = stereoLayer.getStereoSensor();
            
            double axialZ1 = axialSensor1.getGeometry().getPosition().z();
            double stereoZ1 = stereoSensor1.getGeometry().getPosition().z(); 
            
            if (stereoLayer.getLayerNumber() == 1) {
                if (axialSensor1.isTopLayer()) topLayer1Z = axialZ1; 
                else botLayer1Z = axialZ1;
            } else if(stereoLayer.getLayerNumber() == 2) {
                if (axialSensor1.isTopLayer()) topLayer2Z = stereoZ1;
                else botLayer2Z = stereoZ1;
            }
            
            // 2nd stereo pair (layers 3 & 4)
            if (stereoLayer.getLayerNumber() < 3 && stereoLayer.getLayerNumber() > 4) continue;
            
            HpsSiSensor axialSensor2 = stereoLayer.getAxialSensor();
            HpsSiSensor stereoSensor2 = stereoLayer.getStereoSensor();
            
            double axialZ2 = axialSensor2.getGeometry().getPosition().z();
            double stereoZ2 = stereoSensor2.getGeometry().getPosition().z(); 
            
            if (stereoLayer.getLayerNumber() == 3) {
                if (axialSensor2.isTopLayer()) topLayer3Z = axialZ2; 
                else botLayer3Z = axialZ2;
            } else if(stereoLayer.getLayerNumber() == 4) {
                if (axialSensor2.isTopLayer()) topLayer4Z = stereoZ2; 
                else botLayer4Z = stereoZ2;
            }
            
         // 3rd stereo pair (layers 5 & 6)
            if (stereoLayer.getLayerNumber() < 5 && stereoLayer.getLayerNumber() > 6) continue;
            
            HpsSiSensor axialSensor3 = stereoLayer.getAxialSensor();
            HpsSiSensor stereoSensor3 = stereoLayer.getStereoSensor();
            
            double axialZ3 = axialSensor3.getGeometry().getPosition().z();
            double stereoZ3 = stereoSensor3.getGeometry().getPosition().z(); 
            
            if (stereoLayer.getLayerNumber() == 5) {
                if (axialSensor3.isTopLayer()) topLayer5Z = axialZ3; 
                else botLayer5Z = axialZ3;
            } else if(stereoLayer.getLayerNumber() == 6) {
                if (axialSensor3.isTopLayer()) topLayer6Z = stereoZ3; 
                else botLayer6Z = stereoZ3;
            }
            
         // 4th stereo pair (layers 7 & 8)
            if (stereoLayer.getLayerNumber() < 7 && stereoLayer.getLayerNumber() > 8) continue;
            
            HpsSiSensor axialSensor4 = stereoLayer.getAxialSensor();
            HpsSiSensor stereoSensor4 = stereoLayer.getStereoSensor();
            
            double axialZ4 = axialSensor4.getGeometry().getPosition().z();
            double stereoZ4 = stereoSensor4.getGeometry().getPosition().z(); 
            
            if (stereoLayer.getLayerNumber() == 7) {
                if (axialSensor4.isTopLayer()) topLayer7Z = axialZ4; 
                else botLayer7Z = axialZ4;
            } else if(stereoLayer.getLayerNumber() == 8) {
                if (axialSensor4.isTopLayer()) topLayer8Z = stereoZ4; 
                else botLayer8Z = stereoZ4;
            }
            
         // 5th stereo pair (layers 9 & 10)
            if (stereoLayer.getLayerNumber() < 9 && stereoLayer.getLayerNumber() > 10) continue;
            
            HpsSiSensor axialSensor5 = stereoLayer.getAxialSensor();
            HpsSiSensor stereoSensor5 = stereoLayer.getStereoSensor();
            
            double axialZ5 = axialSensor5.getGeometry().getPosition().z();
            double stereoZ5 = stereoSensor5.getGeometry().getPosition().z(); 
            
            if (stereoLayer.getLayerNumber() == 9) {
                if (axialSensor5.isTopLayer()) topLayer9Z = axialZ5; 
                else botLayer9Z = axialZ5;
            } else if(stereoLayer.getLayerNumber() == 10) {
                if (axialSensor5.isTopLayer()) topLayer10Z = stereoZ5; 
                else botLayer10Z = stereoZ5;
            }
            
         // 6th stereo pair (layers 11 & 12)
            if (stereoLayer.getLayerNumber() < 11 && stereoLayer.getLayerNumber() > 12) continue;
            
            HpsSiSensor axialSensor6 = stereoLayer.getAxialSensor();
            HpsSiSensor stereoSensor6 = stereoLayer.getStereoSensor();
            
            double axialZ6 = axialSensor6.getGeometry().getPosition().z();
            double stereoZ6 = stereoSensor6.getGeometry().getPosition().z(); 
            
            if (stereoLayer.getLayerNumber() == 11) {
                if (axialSensor6.isTopLayer()) topLayer11Z = axialZ6; 
                else botLayer11Z = axialZ6;
            } else if(stereoLayer.getLayerNumber() == 12) {
                if (axialSensor6.isTopLayer()) topLayer12Z = stereoZ6; 
                else botLayer12Z = stereoZ6;
            }
            
            
        } // stereoLayers
    } // detectorChanged
    
    @Override
    public void process(EventHeader event) {
    
        // If the event doesn't have the specified collection of tracks, throw
        // an exception.
        if (!event.hasCollection(Track.class, gblTrackCollectionName)) {
            throw new RuntimeException("Track collection " + gblTrackCollectionName + " doesn't exist");
        }
        
        // Get the collection of tracks from the event
        List<Track> tracks = event.get(Track.class, gblTrackCollectionName);
        
        // Loop through all tracks in an event and tweak the track parameters
        for (Track track : tracks) { 
            
           // Get the track state at the target
          // TrackState trackState = track.getTrackStates().get(0);
           
           // *********************************************************************
           // Extrapolate to layers 1-6, then to the ECal, save track states, 
           // and write to the event.
           //
           // *********************************************************************
          
           // get track state at target
           TrackState stateIP = TrackUtils.getTrackStateAtLocation(track, TrackState.AtIP);
           if (stateIP == null) { 
               throw new RuntimeException("IP track state for GBL track was not found");
           }
    
           // TODO Replace AtOther place holders with new trackState locations defined in TrackState class
           
           // Get the track state at L1 (axial), extrapolated from IP
           double layer1Z = stateIP.getTanLambda() > 0 ? topLayer1Z : botLayer1Z; 
           TrackState stateLayer1 = TrackUtils.extrapolateTrackUsingFieldMap(stateIP, extStartPos, layer1Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer1).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer1);

           // Get the track state at L2 (stereo), extrapolated from L1
           double layer2Z = stateLayer1.getTanLambda() > 0 ? topLayer2Z : botLayer2Z; 
           TrackState stateLayer2 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer1, layer1Z, layer2Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer2).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer2);
           
           // Get the track state at L3 (axial)
           double layer3Z = stateLayer2.getTanLambda() > 0 ? topLayer3Z : botLayer3Z; 
           TrackState stateLayer3 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer2, layer2Z, layer3Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer3).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer3);
           
           // Get the track state at L4 (stereo)
           double layer4Z = stateLayer3.getTanLambda() > 0 ? topLayer4Z : botLayer4Z; 
           TrackState stateLayer4 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer3, layer3Z, layer4Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer4).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer4);
           
           // Get the track state at L5 (axial)
           double layer5Z = stateLayer4.getTanLambda() > 0 ? topLayer5Z : botLayer5Z; 
           TrackState stateLayer5 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer4, layer4Z, layer5Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer5).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer5);
           
           // Get the track state at L6 (stereo)
           double layer6Z = stateLayer5.getTanLambda() > 0 ? topLayer6Z : botLayer6Z; 
           TrackState stateLayer6 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer5, layer5Z, layer6Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer6).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer6);
           
           // Get the track state at L7 (axial)
           double layer7Z = stateLayer6.getTanLambda() > 0 ? topLayer7Z : botLayer7Z; 
           TrackState stateLayer7 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer6, layer6Z, layer7Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer7).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer7);
           
           // Get the track state at L8 (stereo)
           double layer8Z = stateLayer7.getTanLambda() > 0 ? topLayer8Z : botLayer8Z; 
           TrackState stateLayer8 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer7, layer7Z, layer8Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer8).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer8);
           
           // Get the track state at L9 (axial)
           double layer9Z = stateLayer8.getTanLambda() > 0 ? topLayer9Z : botLayer9Z; 
           TrackState stateLayer9 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer8, layer8Z, layer9Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer9).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer9);
           
           // Get the track state at L10 (stereo)
           double layer10Z = stateLayer9.getTanLambda() > 0 ? topLayer10Z : botLayer10Z; 
           TrackState stateLayer10 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer9, layer9Z, layer10Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer10).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer10);
         
           // Begin fringe field

           // Get the track state at L11 (axial)
           double layer11Z = stateLayer10.getTanLambda() > 0 ? topLayer11Z : botLayer11Z; 
           TrackState stateLayer11 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer10, layer10Z, layer11Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer11).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer11);
           
           // Get the track state at L12 (stereo)
           double layer12Z = stateLayer11.getTanLambda() > 0 ? topLayer12Z : botLayer12Z; 
           TrackState stateLayer12 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer11, layer11Z, layer12Z, stepSize, bFieldMap);
           ((BaseTrackState) stateLayer12).setLocation(TrackState.AtOther);
           track.getTrackStates().add(stateLayer12);
           
           // Extrapolate to ECal from SVT last SVT layer
           TrackState stateEcalfromL6 = TrackUtils.extrapolateTrackUsingFieldMap(stateLayer12, layer12Z, ecalPosition, stepSize, bFieldMap);
           
           // Replace the existing track state at the Ecal
           int ecalTrackStateIndex = track.getTrackStates().indexOf(TrackUtils.getTrackStateAtLocation(track, TrackState.AtCalorimeter));
           track.getTrackStates().set(ecalTrackStateIndex, stateEcalfromL6);
     
           
        } // loop over tracks
   
        // If the event doesn't have the specified collection of tracks, throw
        // an exception.
        if (!event.hasCollection(Track.class, seedTrackCollectionName)) {
            throw new RuntimeException("Track collection " + seedTrackCollectionName + " doesn't exist");
        }
        
        // Get the collection of seed tracks from the event and force 
        // the recomputation of the momentum.  This is a hack to force
        // the persistence of the momentum, otherwise, a bogus momentum
        // value is used.
        List<Track> seedTracks = event.get(Track.class, seedTrackCollectionName);
        for (Track seedTrack : seedTracks) { 
            
            // Get the track state at the target
            TrackState trackState = seedTrack.getTrackStates().get(0);
            
            // Force re-computation of momentum using the correct B-field, 
            // otherwise, a bogus value is returned. 
            ((BaseTrack) seedTrack).setTrackParameters(seedTrack.getTrackStates().get(0).getParameters(), bField);
        }
        
    }
}
