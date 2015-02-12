package org.hps.recon.ecal;

import java.util.Map;

import org.hps.conditions.database.TableConstants;
import org.hps.conditions.ecal.EcalChannel;
import org.hps.conditions.ecal.EcalChannelConstants;
import org.hps.conditions.ecal.EcalConditions;
import org.lcsim.conditions.ConditionsManager;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.RawCalorimeterHit;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.base.BaseRawCalorimeterHit;
import org.lcsim.geometry.Detector;

/**
 * This class is used to convert {@link org.lcsim.event.RawCalorimeterHit}
 * objects to {@link org.lcsim.event.CalorimeterHit} objects with energy
 * information. It has methods to convert pedestal subtracted ADC counts to
 * energy.
 *
 * @author Sho Uemura <meeg@slac.stanford.edu>
 * @author Jeremy McCormick <jeremym@slac.stanford.edu>
 * @author Andrea Celentano <andrea.celentano@ge.infn.it>
 */
public class EcalRawConverter {

    private boolean useTimeWalkCorrection = false;
    private boolean useRunningPedestal = false;
    private boolean constantGain = false;
    private double gain;
    private boolean use2014Gain = true;

    private EcalConditions ecalConditions = null;

    public EcalRawConverter() {
    }

    public void setGain(double gain) {
        constantGain = true;
        this.gain = gain;
    }

    public void setUse2014Gain(boolean use2014Gain) {
        this.use2014Gain = use2014Gain;
    }

    public void setUseRunningPedestal(boolean useRunningPedestal) {
        this.useRunningPedestal=useRunningPedestal;
    }
    
    public void setUseTimeWalkCorrection(boolean useTimeWalkCorrection) {
        this.useTimeWalkCorrection=useTimeWalkCorrection;
    }
  
    /*
     * NAB 2015/02/11 
     * Choose whether to use static pedestal from database or running pedestal.
     * This can only used for Mode-7 data.
     */
    public double getMode7Pedestal(EventHeader event,RawCalorimeterHit hit)
    {
        if (useRunningPedestal) {
            if (event.hasItem("EcalRunningPedestals")) {
                Map<EcalChannel, Double> runningPedMap=
                        (Map<EcalChannel, Double>)
                        event.get("EcalRunningPedestals");
                EcalChannel chan = ecalConditions.getChannelCollection().
                        findGeometric(hit.getCellID());
                //System.err.println(" %%%%%%%%%%%%%%%%% "+chan.getChannelId()+" "+runningPedMap.get(chan));
                if (!runningPedMap.containsKey(chan)){
                    System.err.println("************** Missing Pedestal");
                } else {
                    return runningPedMap.get(chan);
                }
            } else {
                System.err.println("*****************************************************************");
                System.err.println("**  You Requested a Running Pedestal, but it is NOT available. **");
                System.err.println("**     Reverting to the database. Only printing this ONCE.     **");
                System.err.println("*****************************************************************");
                useRunningPedestal = false;
            }
        }
        return findChannel(hit.getCellID()).getCalibration().getPedestal();
    }
    public short sumADC(RawTrackerHit hit) {
        EcalChannelConstants channelData = findChannel(hit.getCellID());
        double pedestal = channelData.getCalibration().getPedestal();
        short sum = 0;
        short samples[] = hit.getADCValues();
        for (int isample = 0; isample < samples.length; ++isample) {
            sum += (samples[isample] - pedestal);
        }
        return sum;
    }

    public CalorimeterHit HitDtoA(RawTrackerHit hit) {
        double time = hit.getTime();
        long id = hit.getCellID();
        double rawEnergy = adcToEnergy(sumADC(hit), id);
        return CalorimeterHitUtilities.create(rawEnergy, time, id);
    }

    /*
     * This HitDtoA is for Mode-3 data at least, but definitely not Mode-7.
     * A time-walk correction can be applied.  (NAB 2015/02/11).
     */
    public CalorimeterHit HitDtoA(RawCalorimeterHit hit, int window, double timeOffset) {
        if (hit.getTimeStamp() % 64 != 0) {
            System.out.println("unexpected timestamp " + hit.getTimeStamp());
        }
        double time = hit.getTimeStamp() / 16.0;
        long id = hit.getCellID();
        EcalChannelConstants channelData = findChannel(id);
        double adcSum = hit.getAmplitude() - window * channelData.getCalibration().getPedestal();
        double rawEnergy = adcToEnergy(adcSum, id);
        if (useTimeWalkCorrection) {
           time = EcalTimeWalk.correctTimeWalk(time,rawEnergy);
        }
        return CalorimeterHitUtilities.create(rawEnergy, time + timeOffset, id);
    }

    /*
     * This HitDtoA is exclusively for Mode-7 data, hence the GenericObject parameter.
     * The decision to call this method is made in EcalRawConverterDriver based on the
     * format of the input EVIO data.  EventHeader is also passed in order to allow access
     * to running pedestals, which is only applicable to Mode-7 data.  (NAB, 2015/02/11)
     */
    public CalorimeterHit HitDtoA(EventHeader event,RawCalorimeterHit hit, GenericObject mode7Data, int window, double timeOffset) {
        double time = hit.getTimeStamp() / 16.0; //timestamps use the full 62.5 ps resolution
        long id = hit.getCellID();
        double adcSum = hit.getAmplitude() - window * getMode7Pedestal(event,hit);
        double rawEnergy = adcToEnergy(adcSum, id);        
        return CalorimeterHitUtilities.create(rawEnergy, time + timeOffset, id);
    }

    public RawCalorimeterHit HitAtoD(CalorimeterHit hit, int window) {
        int time = (int) (Math.round(hit.getTime() / 4.0) * 64.0);
        long id = hit.getCellID();
        // Get the channel data.
        EcalChannelConstants channelData = findChannel(id);
        int amplitude;
        if (constantGain) {
            amplitude = (int) Math.round((hit.getRawEnergy() / ECalUtils.MeV) / gain + window * channelData.getCalibration().getPedestal());
        } else {
            amplitude = (int) Math.round((hit.getRawEnergy() / ECalUtils.MeV) / channelData.getGain().getGain() + window * channelData.getCalibration().getPedestal());
        }
        RawCalorimeterHit h = new BaseRawCalorimeterHit(id, amplitude, time);
        return h;
    }

    /*
     * return energy (units of GeV) corresponding to the ADC sum and crystal ID
     */
    private double adcToEnergy(double adcSum, long cellID) {

        // Get the channel data.
        EcalChannelConstants channelData = findChannel(cellID);

        if (use2014Gain) {
            if (constantGain) {
                return adcSum * ECalUtils.gainFactor * ECalUtils.ecalReadoutPeriod;
            } else {
                return channelData.getGain().getGain() * adcSum * ECalUtils.gainFactor * ECalUtils.ecalReadoutPeriod; // should not be used for the moment (2014/02)
            }
        } else {
            if (constantGain) {
                return gain * adcSum * ECalUtils.MeV;
            } else {
                return channelData.getGain().getGain() * adcSum * ECalUtils.MeV; //gain is defined as MeV/integrated ADC
            }
        }
    }

    /**
     * Must be set when an object EcalRawConverter is created.
     *
     * @param detector (long)
     */
    public void setDetector(Detector detector) {
        // ECAL combined conditions object.
        ecalConditions = ConditionsManager.defaultInstance()
                .getCachedConditions(EcalConditions.class, TableConstants.ECAL_CONDITIONS).getCachedData();
    }

    /**
     * Convert physical ID to gain value.
     *
     * @param cellID (long)
     * @return channel constants (EcalChannelConstants)
     */
    public EcalChannelConstants findChannel(long cellID) {
        return ecalConditions.getChannelConstants(ecalConditions.getChannelCollection().findGeometric(cellID));
    }    
}
