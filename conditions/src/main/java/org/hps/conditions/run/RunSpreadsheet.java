package org.hps.conditions.run;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

/**
 * A simple representation of the 2015 run spreadsheet (runs from 3/28 to 5/19) read from an exported CSV file as a list of records.
 * <p>
 * Master copy of the spreadsheet is located at 
 * <a href="https://docs.google.com/spreadsheets/d/1l1NurPpsmpgZKgr1qoQpLQBBLz1sszLz4xZF-So4xs8/edit#gid=43855609">HPS_Runs_2015</a>.
 * <p>
 * The rows are accessible as raw CSV data through the Apache Commons CSV library, and this data must be manually cleaned up and converted 
 * to the correct data type before being inserted into the conditions database.
 *
 * @author Jeremy McCormick
 */
public final class RunSpreadsheet {

    /**
     * The column headers.
     */
    private static String[] HEADERS = {"run", "date", "start_time", "end_time", "to_tape", "n_events", "trigger_rate", "target", "beam_current",
        "beam_x", "beam_y", "trigger_config", "ecal_fadc_mode", "ecal_fadc_thresh", "ecal_fadc_window", "ecal_cluster_thresh_seed", "ecal_cluster_thresh_cluster",
        "ecal_cluster_window_hits", "ecal_cluster_window_pairs", "ecal_scalers_fadc", "ecal_scalers_dsc", "svt_y_position", "svt_offset_phase", "svt_offset_time",
        "ecal_temp", "ecal_lv_current", "notes"};

    /**
     * Read the CSV file from the command line and print the data to the terminal (just a basic test).
     *
     * @param args the command line arguments
     */
    public static void main(final String args[]) throws Exception {
        final RunSpreadsheet runSpreadsheet = new RunSpreadsheet(new File(args[0]));
        for (final CSVRecord record : runSpreadsheet.getRecords()) {
            try {
                System.out.print("start date: " + parseStartDate(record) + ", ");
                System.out.print("end date: " + parseEndDate(record) + ", ");
                System.out.print(record);
                System.out.println();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * The input CSV file.
     */
    private final File file;

    /**
     * The list of records read from the CSV file.
     */
    private List<CSVRecord> records;

    /**
     * Read the data from a CSV export of the run spreadsheet.
     *
     * @param file the CSV file
     */
    public RunSpreadsheet(final File file) {
        this.file = file;
        try {
            this.fromCsv(this.file);
        } catch (final Exception e) {
            throw new RuntimeException();
        }
    }

    /**
     * Find the record for a run in the list of records.
     *
     * @param run the run number
     * @return the <code>CSVRecord</code> or <code>null</code> if not found
     */
    public CSVRecord findRun(final int run) {
        for (final CSVRecord record : records) {
            try {
                if (run == Integer.parseInt(record.get("run_number"))) {
                    return record;
                }
            } catch (final NumberFormatException e) {
                e.printStackTrace();
            }
        }
        return null;
    }

    /**
     * Read all the records from the CSV file.
     *
     * @param file the CSV file
     * @throws FileNotFoundException if file does not exist
     * @throws IOException if there is an IO error
     */
    private void fromCsv(final File file) throws FileNotFoundException, IOException {

        final FileReader reader = new FileReader(file);
        final CSVFormat format = CSVFormat.DEFAULT.withHeader(HEADERS);

        final CSVParser parser = new CSVParser(reader, format);

        records = parser.getRecords();

        // Remove first three rows of headers.
        records.remove(0);
        records.remove(0);
        records.remove(0);

        parser.close();
    }

    /**
     * Get the list of records read from the CSV file.
     *
     * @return the list of records read from the CSV file
     */
    public List<CSVRecord> getRecords() {
        return records;
    }
    
    private static final SimpleDateFormat DATE_FORMAT = new SimpleDateFormat("MM/dd/yyyy H:mm"); 
    
    private static Date parseStartDate(CSVRecord record) throws ParseException {
        return DATE_FORMAT.parse(record.get("date") + " " + record.get("start_time"));
    }
    
    private static Date parseEndDate(CSVRecord record) throws ParseException {
        return DATE_FORMAT.parse(record.get("date") + " " + record.get("end_time"));
    }
    
    private static int parseRunNumber(CSVRecord record) throws NumberFormatException {
        return Integer.parseInt(record.get("run"));
    }
    
    public static class RunData {
        
        private int run;
        private Date startDate;
        private Date endDate;
        private CSVRecord record;
        
        RunData(CSVRecord record) throws NumberFormatException {
            this.record = record;
            run = parseRunNumber(this.record);
            try {
                startDate = RunSpreadsheet.parseStartDate(this.record);
            } catch (ParseException e) {                
            }
            try {
                endDate = RunSpreadsheet.parseEndDate(this.record);
            } catch (ParseException e) {                
            }
        }
        
        public int getRun() {
            return run;
        }
        
        public Date getStartDate() {
            return startDate;
        }
        
        public Date getEndDate() {
            return endDate;
        }      
        
        public String toString() {
            return "RunData { run: " + run + ", startDate: " + startDate + ", endDate: " + endDate + " }";
        }
        
        public CSVRecord getRecord() {
            return record;
        }
    }
    
    @SuppressWarnings("serial")
    public static class RunMap extends LinkedHashMap<Integer, RunData> {
        
        private void addRunData(RunData runData) {
            this.put(runData.getRun(), runData);
        }
    }
    
    public RunMap getRunMap() {
        RunMap runMap = new RunMap();
        for (final CSVRecord record : getRecords()) {
            try {
                runMap.addRunData(new RunData(record));
            } catch (NumberFormatException e) {
                e.printStackTrace();
            }
        }
        return runMap;
    }    
}