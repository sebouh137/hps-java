package org.hps.record.evio.crawler;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import org.hps.record.evio.EvioEventProcessor;
import org.jlab.coda.jevio.EvioEvent;
import org.jlab.coda.jevio.EvioException;
import org.jlab.coda.jevio.EvioReader;
import org.lcsim.util.log.LogUtil;

public final class RunProcessor {

    private static final Logger LOGGER = LogUtil.create(RunProcessor.class);

    List<EvioEventProcessor> processors = new ArrayList<EvioEventProcessor>();

    RunSummary runSummary;

    RunProcessor(final RunSummary runSummary) {
        this.runSummary = runSummary;
    }

    void addProcessor(final EvioEventProcessor processor) {
        this.processors.add(processor);
        LOGGER.config("added processor: " + processor.getClass().getSimpleName());
    }

    List<EvioEventProcessor> getProcessors() {
        return this.processors;
    }

    void process() throws Exception {
        
        // Run the start of job hooks.
        for (final EvioEventProcessor processor : this.processors) {
            processor.startJob();
        }
         
        
        
        // Process the files.
        for (final File file : this.runSummary.getFiles()) {
            process(file);
        }
        
        // Run the end of job hooks.
        for (final EvioEventProcessor processor : this.processors) {
            processor.endJob();
        }
    }

    private void process(final File file) throws EvioException, IOException, Exception {
        EvioReader reader = null;
        try {
            // Open with wrapper method which will use the cached file path if necessary.
            reader = EvioFileUtilities.open(file);
            
            // Compute event count for the file and store the value.
            this.runSummary.getFiles().computeEventCount(reader, file);
            
            // Process the events using the list of EVIO processors.
            EvioEvent event = null;
            while ((event = reader.parseNextEvent()) != null) {
                for (final EvioEventProcessor processor : this.processors) {
                    processor.process(event);
                }
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }
   
}
