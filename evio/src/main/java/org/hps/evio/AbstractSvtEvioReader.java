package org.hps.evio;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.hps.record.svt.SvtEvioUtils;
import org.hps.record.svt.SvtHeaderDataInfo;
import org.hps.record.svt.SvtEvioExceptions.SvtEvioHeaderException;
import org.hps.record.svt.SvtEvioExceptions.SvtEvioReaderException;
import org.hps.util.Pair;
import org.jlab.coda.jevio.BaseStructure;
import org.jlab.coda.jevio.EvioEvent;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.base.BaseRawTrackerHit;
import org.lcsim.geometry.Subdetector;
import org.lcsim.lcio.LCIOUtil;

/**
 * Abstract SVT EVIO reader used to convert SVT bank sample blocks to
 * {@link RawTrackerHit}s.
 * 
 * @author Omar Moreno <omoreno1@ucsc.edu>
 * @author Per Hansson Adrian <phansson@slac.stanford.edu>
 *
 */
public abstract class AbstractSvtEvioReader extends EvioReader {
    
    public static final String SVT_HEADER_COLLECTION_NAME = "SvtHeaders";
    
    // Initialize the logger
    public static final Logger LOGGER = Logger.getLogger(AbstractSvtEvioReader.class.getPackage().getName());
    
    // A Map from DAQ pair (FPGA/Hybrid or FEB ID/FEB Hybrid ID) to the
    // corresponding sensor
    protected Map<Pair<Integer /* FPGA */, Integer /* Hybrid */>,
                  HpsSiSensor /* Sensor */> daqPairToSensor 
                      = new HashMap<Pair<Integer, Integer>, HpsSiSensor>();
    
    // Flag indicating whether the DAQ map has been setup
    protected boolean isDaqMapSetup = false;

    // Collections and names
    private static final String SVT_HIT_COLLECTION_NAME = "SVTRawTrackerHits";

    // Constants
    private static final String SUBDETECTOR_NAME = "Tracker";
    private static final String READOUT_NAME = "TrackerHits";

    /**
     *  Get the minimum SVT ROC bank tag in the event.
     *
     *  @return Minimum SVT ROC bank tag
     */
    abstract protected int getMinRocBankTag(); 
    
    /**
     *  Get the maximum SVT ROC bank tag in the event.
     *
     *  @return Maximum SVT ROC bank tag
     */
    abstract protected int getMaxRocBankTag(); 
    
    /**
     *  Get the minimum SVT ROC bank tag in the event.
     *
     *  @return Minimum SVT ROC bank tag
     */
    abstract protected int getMinDataBankTag(); 
    
    /**
     *  Get the maximum SVT ROC bank tag in the event.
     *
     *  @return Maximum SVT ROC bank tag
     */
    abstract protected int getMaxDataBankTag(); 

    /**
     *  Get the SVT ROC bank number of the bank encapsulating the SVT samples.
     * 
     *  @return SVT ROC bank number 
     */
    abstract protected int getRocBankNumber(); 
    
    /**
     *  Get the number of 32 bit integers composing the data block header
     *
     *  @return The header length
     */
    abstract protected int getDataHeaderLength();

    /**
     *  Get the number of 32 bit integers composing the data block tail (the 
     *  data inserted after all sample blocks in a data block)
     * 
     *  @return The tail length 
     */
    abstract protected int getDataTailLength();

    /**
     *  A method to setup a mapping between a DAQ pair 
     *  (FPGA/Hybrid or FEB ID/FEB Hybrid ID) and the corresponding sensor.
     *
     *  @param subdetector - The tracker {@link Subdetector} object
     */
    // TODO: This can probably be done when the conditions are loaded.
    abstract protected void setupDaqMap(Subdetector subdetector);

    /**
     *  Get the sensor associated with a set of samples  
     *
     *  @param data - sample block of data
     *  @return The sensor associated with a set of sample 
     */
    abstract protected HpsSiSensor getSensor(int[] data);

    /**
     * Check whether the samples are valid
     * 
     * @param data : sample block of data
     * @return true if the samples are valid, false otherwise
     */
    abstract protected boolean isValidSampleSet(int[] data);
    
    /**
     *  Process an EVIO event and extract all information relevant to the SVT.
     *  
     *  @param event - EVIO event to process
     *  @param lcsimEvent - LCSim event to put collections into 
     *  @return true if the EVIO was processed successfully, false otherwise 
     * @throws SvtEvioReaderException 
     */
    public boolean processEvent(EvioEvent event, EventHeader lcsimEvent) throws SvtEvioReaderException {
        return this.makeHits(event, lcsimEvent);
    }

    /**
     *  Make {@link RawTrackerHit}s out of all sample sets in an SVT EVIO bank
     *  and put them into an LCSim event.
     *
     *  
     *  @param event - EVIO event to process
     *  @param lcsimEvent - LCSim event to put collections into 
     *  @return true if the raw hits were created successfully, false otherwise 
     * @throws SvtEvioReaderException 
     */
    @Override
    public boolean makeHits(EvioEvent event, EventHeader lcsimEvent) throws SvtEvioReaderException {

        LOGGER.finest("Physics Event: " + event.toString());
        
        // Retrieve the data banks encapsulated by the physics bank.  The ROC
        // bank range is set in the subclass.
        List<BaseStructure> dataBanks = SvtEvioUtils.getDataBanks(event, this.getMinRocBankTag(), this.getMaxRocBankTag(), this.getMinDataBankTag(), this.getMaxDataBankTag());
        
        // Return false if data banks weren't found
        if (dataBanks.isEmpty()) return false;  
    
        // Setup the DAQ map if it's not setup
        if (!this.isDaqMapSetup)
            this.setupDaqMap(lcsimEvent.getDetector().getSubdetector(
                    SUBDETECTOR_NAME));

        List<RawTrackerHit> rawHits = new ArrayList<RawTrackerHit>();
        List<SvtHeaderDataInfo> headers = new ArrayList<SvtHeaderDataInfo>();

        LOGGER.finest("Total data banks found: " + dataBanks.size());

            // Loop over all of the data banks contained by the ROC banks and 
        // processed them
        for (BaseStructure dataBank : dataBanks) {

            LOGGER.finest("Processing data bank: " + dataBank.toString());

            // Get the int data encapsulated by the data bank
            int[] data = dataBank.getIntData();
            LOGGER.finest("Total number of integers contained by the data bank: " + data.length);

            // Check that a complete set of samples exist
            int sampleCount = data.length - this.getDataHeaderLength()
                    - this.getDataTailLength();
            LOGGER.finest("Total number of  samples: " + sampleCount);
            if (sampleCount % 4 != 0) {
                throw new SvtEvioReaderException("[ "
                        + this.getClass().getSimpleName()
                        + " ]: Size of samples array is not divisible by 4");
            }

            // extract header and tail information
            SvtHeaderDataInfo headerData = this.extractSvtHeader(dataBank.getHeader().getNumber(), data);

            // Check that the multisample count is consistent
            this.checkSvtSampleCount(sampleCount, headerData);

            // Add header to list
            headers.add(headerData);

            // Store the multisample headers
            // Note that the length is not known but can't be longer than the multisample count
            // in other words the data can be only header multisamples for example.
            int multisampleHeaderData[] = new int[sampleCount];
            int multisampleHeaderIndex = 0;

            LOGGER.finest("sampleCount " + sampleCount);

            List<int[]> multisampleList = SvtEvioUtils.getMultisamples(data, sampleCount, this.getDataHeaderLength());
            // Loop through all of the samples and make hits
            for (int[] samples:multisampleList) {
                if (SvtEvioUtils.isMultisampleHeader(samples)) {
                    LOGGER.finest("this is a header multisample for apv " + SvtEvioUtils.getApvFromMultiSample(samples) + " ch " + SvtEvioUtils.getChannelNumber(samples));
                    // Extract data words from multisample header and update index
                    multisampleHeaderIndex += this.extractMultisampleHeaderData(samples, multisampleHeaderIndex, multisampleHeaderData);
                } else {
                    LOGGER.finest("this is a data multisample for apv " + SvtEvioUtils.getApvFromMultiSample(samples) + " ch " + SvtEvioUtils.getChannelNumber(samples));
                }

                // If a set of samples is associated with an APV header or tail, skip it
                if (!this.isValidSampleSet(samples)) {
                    continue;
                }
                rawHits.add(this.makeHit(samples));
            }

            LOGGER.finest("got " + multisampleHeaderIndex + " multisampleHeaderIndex for " + sampleCount + " sampleCount");

            // add multisample header tails to header data object
            this.setMultiSampleHeaders(headerData, multisampleHeaderIndex, multisampleHeaderData);

        }
        
        LOGGER.finest("Total number of RawTrackerHits created: " + rawHits.size());

        // Turn on 64-bit cell ID.
        int flag = LCIOUtil.bitSet(0, 31, true);
        // Add the collection of raw hits to the LCSim event
        lcsimEvent.put(SVT_HIT_COLLECTION_NAME, rawHits, RawTrackerHit.class, flag, READOUT_NAME);

        // Process SVT headers
        this.processSvtHeaders(headers, lcsimEvent);
        
        return true;
    }

    
    
    /**
     * Process the headers that were extracted from the SVT data. 
     * @param headers - list of all headers
     * @param lcsimEvent - the current LCSIM event being processed
     * @throws SvtEvioHeaderException
     */
    protected abstract void processSvtHeaders(List<SvtHeaderDataInfo> headers, EventHeader lcsimEvent) throws SvtEvioHeaderException;
    
    /**
     * Extract the header information and store it in a {@link SvtHeaderDataInfo} object.
     * @param num - bank num (ROC id)
     * @param data - SVT data block.
     * @return the {@link SvtHeaderDataInfo}.
     */
    protected SvtHeaderDataInfo extractSvtHeader(int num, int[] data) {
        // Extract the header information
        int svtHeader = SvtEvioUtils.getSvtHeader(data);
        // Extract the tail information
        int svtTail = SvtEvioUtils.getSvtTail(data);
        return new SvtHeaderDataInfo(num, svtHeader, svtTail);

    }
    
    /**
     * Copy the multisample header data for the samples into a long array.
     * @param samples
     * @param index
     * @param multisampleHeaderData
     * @return the length of the extracted samples or 0 if not a multisample header
     */
    protected int extractMultisampleHeaderData(int[] samples, int index, int[] multisampleHeaderData) {
        LOGGER.finest("extractMultisampleHeaderData: index " + index);
        if( SvtEvioUtils.isMultisampleHeader(samples) && !SvtEvioUtils.isMultisampleTail(samples) ) {
            LOGGER.finest("extractMultisampleHeaderData: this is a multisample header so add the words to index " + index);
            System.arraycopy(samples, 0, multisampleHeaderData, index, samples.length);
            return samples.length;
        } else {
            LOGGER.finest("extractMultisampleHeaderData: this is a NOT multisample header ");
            return 0;
        }
    }
    
    /**
     * Checks that the SVT header data count is consistent with the bank size.
     * @param sampleCount - sample count from the size.
     * @param headerData - header extracted from the bank.
     * @throws SvtEvioHeaderException
     */
    protected void checkSvtSampleCount(int sampleCount, SvtHeaderDataInfo headerData) throws SvtEvioHeaderException {
        if( sampleCount != SvtEvioUtils.getSvtTailMultisampleCount(headerData.getTail())*4)
            throw new SvtEvioHeaderException("ROC " + headerData.getNum() + " multisample count " + sampleCount + " is not consistent with bank size " + SvtEvioUtils.getSvtTailMultisampleCount(headerData.getTail()));
    }
    
    /**
     * Add the multisample headers to the {@link SvtHeaderDataInfo} object.
     * @param headerData - object to add multisample headers to.
     * @param n - number of multisample headers
     * @param multisampleHeaders - multisample headers to copy
     */
    protected void setMultiSampleHeaders(SvtHeaderDataInfo headerData, int n, int[] multisampleHeaders) {
        //copy out the headers that are non-zero
        int[] vals = new int[n];
        System.arraycopy(multisampleHeaders, 0, vals, 0, n);
        //logger.info("setMultiSampleHeaders: adding " + vals.length + " multisample headers");
        headerData.setMultisampleHeaders(vals);
    }
    
       
    
    
    /**
     *  Make a {@link RawTrackerHit} from a set of samples.
     * 
     *  @param data : sample block of data
     *  @return A raw hit
     */
    protected abstract RawTrackerHit makeHit(int[] data); 
    
    /**
     *  Make a {@link RawTrackerHit} from a set of samples.
     * 
     *  @param data : Sample block of data
     *  @param channel : Channel number associated with these samples
     *  @return A raw hit
     */
    protected RawTrackerHit makeHit(int[] data, int channel) { 

        // Get the sensor associated with this sample
        HpsSiSensor sensor = this.getSensor(data);
        //logger.fine(sensor.toString());
        
        // Use the channel number to create the cell ID
        long cellID = sensor.makeChannelID(channel);
        
        // Set the hit time.  For now this will be zero
        int hitTime = 0;
    
        // Create and return a RawTrackerHit
        return new BaseRawTrackerHit(hitTime, cellID, SvtEvioUtils.getSamples(data), null, sensor);
    }
}
