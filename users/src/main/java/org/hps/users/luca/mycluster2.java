/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.hps.users.luca;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.util.Driver;

/**
 * 
 * @author Luca Colaneri 
 */
public class mycluster2 extends Driver {
    int posx=0;
    int posy=6;
    int Clustercount=0;
    protected String clusterCollectionName = "EcalClusters";
    
    ArrayList<Cluster> ClusterQueue;
    List<Cluster> goodCluster;
 
    private FileWriter writer;
    //private FileWriter writer2;
    String outputFileName = "ClusterInfoTutti.txt";
   // String outputFileName2 = "ClusterEnePos2.txt";

 // AIDA aida = AIDA.defaultInstance();

   
 
   
 //   IHistogram2D MCendpoint = aida.histogram2D("endpoint MCParticles ",700,-350,350,200,-100,100);
  // IHistogram1D ClusterEne = aida.histogram1D("energia cluster", 100, 0.0,3.0);
   
   
   
   double[] position;
    public void setClusterPosition (int posix, int posiy){
    this.posx=posix;
    this.posy=posiy;
    }
   
    public void setClusterCollectionName(String clusterCollectionName) {
        this.clusterCollectionName = clusterCollectionName;
    }
 
   public void setOutputFileName(String outputFileName){
this.outputFileName = outputFileName;
}
  /*
   *
   *
   *
   */
   
   @Override   
public void startOfData(){

 
    ClusterQueue=new ArrayList<Cluster>();
    goodCluster = new ArrayList<Cluster>();
    
    
    try{
    //initialize the writers
    writer=new FileWriter(outputFileName,true);
 //   writer2=new FileWriter(outputFileName2);
    //Clear the files
    //writer.write("");
  //  writer2.write("");
}
catch(IOException e ){
System.err.println("Error initializing output file for event display.");
}
}
  
@Override
public void endOfData(){
System.out.println("Ho contato" + Clustercount + "clusters \n");
    ClusterAnalyzer();
    try{
//close the file writer.
    writer.close();
   // writer2.close();
    }
catch(IOException e){
    System.err.println("Error closing utput file for event display.");
}
} 
   
 @Override  
 public void process (EventHeader event){
   
          
           
           
     
     
     if(event.hasCollection(Cluster.class, "EcalClusters")) {
            
         fillQueue(event.get(Cluster.class, clusterCollectionName));
          
            
            }       
    
}
 /**
 * Puts cluster collected from events in the ClusterQueue
 */
 public void fillQueue (List<Cluster> ecalClusters){
 for(Cluster cluster : ecalClusters)
 {
     
     ClusterQueue.add(cluster);
 }
 
 
 }
 /**
  * For each crystal, looks for clusters that hit that clystar, if it is an isolated cluster, it's put in goodclusterqueue
  */
 public void ClusterAnalyzer(){
 

 ///cerca i cluster nella posizione che ci interessa poi chiama la funzione che decide se sono "isolati"
   //System.out.println("Sta partendo il for sulla Queue \n");
 for(int y=-5;y<6;y++){
     for(int x=-23;x<24;x++){
      posx=x;
      posy=y;
         
         
    for(int i=0; i< ClusterQueue.size(); i++){
    if(ClusterQueue.get(i).getCalorimeterHits().get(0).getIdentifierFieldValue("ix")== posx && ClusterQueue.get(i).getCalorimeterHits().get(0).getIdentifierFieldValue("iy")==posy && ClusterQueue.get(i).getEnergy()>1.6){
        // System.out.println("sono nel primo if, ovvero ci sono cluster dove voglio io");
           if(ClusterChecker(i)){
      goodCluster.add(ClusterQueue.get(i));}
      }
     
     
    }
 
 
 
 }
 }
 
 
 int id;
 for(Cluster cluster :goodCluster){
   Clustercount++;
   id=getCrystal(cluster);
     try{
     writer.append(id + " " + cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix") + " " + cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy")+ " " + cluster.getEnergy()+ " " + cluster.getSize() + " " + cluster.getCalorimeterHits().get(0).getRawEnergy() + " ");
     /*for(CalorimeterHit hit : cluster.getCalorimeterHits())
     {writer.append(hit.getRawEnergy()+ " ");
       }*/
     writer.append("\n");
     }
     
   catch(IOException e ){System.err.println("Error writing tooutput for event display");}   
     
 }
 // System.out.println("Ci sono " + Clustercount + "cluster! \n" );
 }
 /**
  * Check if the cluster is isolaterd checking if there are clusters near it in time and in space
  * @param pos
  * @return 
  */
 
public boolean ClusterChecker (int pos){
//System.out.println("Sono nel clustercheck! \n");
    
boolean check=true;
    double deltaT=0;
    for (int i=(pos-5);i<pos+5;i++){
        ///controlla che non esca da array
        if(i>=0 && i < ClusterQueue.size()){
            deltaT=Math.abs(ClusterQueue.get(i).getCalorimeterHits().get(0).getTime()-ClusterQueue.get(pos).getCalorimeterHits().get(0).getTime());
    ///controlla la differenza temporale che deve essere >80ns, se è maggiore, check=true, se è minore, checka la posizione
     if(deltaT<80)
     {check=PositionChecker(i);}
            
        
        }
    }
return check;

}
      
 public boolean PositionChecker(int i){
       boolean test=true;

//controllo su posizione orizzontale
     if(ClusterQueue.get(i).getCalorimeterHits().get(0).getIdentifierFieldValue("ix")< posx-2 || ClusterQueue.get(i).getCalorimeterHits().get(0).getIdentifierFieldValue("ix")> posx+2){
        //se orizzonale soddisfatto, controllo verticale
         if(ClusterQueue.get(i).getCalorimeterHits().get(0).getIdentifierFieldValue("iy")< posy-2 || ClusterQueue.get(i).getCalorimeterHits().get(0).getIdentifierFieldValue("ix")> posy+2){
         test= true;
         }
     
     }
     
     else{
     test= false;}
     
     return test;
 }
 
 public int getCrystal (Cluster cluster){
 int x,y,id=0;
 x= (-1)*cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
 y= cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");
 
 if(y==5){
 if(x<0)
 {id=x+24;}
 else id= x+23;
 }
 
 else if(y==4)
 {if(x<0){
  id=x+70;}
 else id=x+69;}
 
 else if(y==3)
 {if(x<0){
  id=x+116;}
 else id=x+115;}
 
 else if(y==2)
 {if(x<0){
  id=x+162;}
 else id=x+161;}
 
 else if(y==1)
 {x=-x;
     if(x>0){
  id=-x+208;}
 else if(x==-1){id=208;}
 else if(x<-1) id=-x+198;}
 
  else if(y==-1)
 {x=-x;
     if(x>0){
  id=-x+245;}
 else if(x==-1 )id=245;
 else if(x<-1){id=-x+257;}}
 
 else if(y==-2)
 {if(x<0){
  id=x+282;}
 else id=x+281;}
 
  else if(y==-3)
 {if(x<0){
  id=x+328;}
 else id=x+327;}
 
 else if(y==-4)
 {if(x<0){
  id=x+374;}
 else id=x+373;}
 
 else if(y==-5)
 {if(x<0){
  id=x+420;}
 else id=x+419;}
 
 return id;
 
 }
 
 } //chiusura driver  
    
    
    
