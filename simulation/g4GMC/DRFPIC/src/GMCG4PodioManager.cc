#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "GMCG4PodioManager.hh"
#include "G4AutoLock.hh"

/// We don't do multithreading at the moment in the full idea simulation. Still, I am leaving the code in for future developments

/*#ifdef G4MULTITHREADED
#include "G4Threading.hh"
#endif*/

#include <string>

GMCG4PodioManager * GMCG4PodioManager::m_inst_ = NULL;  

GMCG4PodioManager::GMCG4PodioManager():
  m_filename_prefix("simhits."),
  m_filename_suffix("podio.root")
{}

GMCG4PodioManager * GMCG4PodioManager::Instance()
{
  if (m_inst_ == NULL){
    m_inst_ = new GMCG4PodioManager();
  }
  return m_inst_;
}

podio::EventStore * GMCG4PodioManager::GetEvtStore()
{
  int threadId =  -10; // some default value different from -1, which is reserved for master in multithread mode

  // Not doing multithreading at the moment. 
  
  /*#ifdef G4MULTITHREADED
  threadId = G4Threading::G4GetThreadId();
  #endif*/

  if (threadId == -1) {// this is the master, nothing to be done
    return NULL;
  }

  // Look if this thread has already a store in the map. If not, create it
  if (m_map_store.find(threadId)  == m_map_store.end()){

    //    G4AutoLock lock(&GMCG4ActionMutex);

    // If it is not there, we need to create it

    podio::EventStore * l_evtstore = new podio::EventStore();

    // build the file name

    G4String filename_id = "";

    if (threadId != -10){
      filename_id = "_";
      filename_id += std::to_string(threadId);
      filename_id += "_";
    }

    G4String filename = m_filename_prefix + filename_id + m_filename_suffix;

    podio::ROOTWriter * l_writer = new podio::ROOTWriter(filename,l_evtstore);

    G4cout << "Podio output file name " << filename << " generated by threadId " << threadId << G4endl;
    
    m_map_store[threadId] = l_evtstore;
    m_map_writer[threadId] = l_writer;
  }

  return m_map_store[threadId];
}

podio::ROOTWriter * GMCG4PodioManager::GetWriter()
{

  podio::EventStore * l_eventstore = this->GetEvtStore();
  if (l_eventstore == NULL) return NULL;
  int threadId =  -10; // some default value different from -1, which is reserved for master in multithread mode

  // we are not doing multithreading
  
  /*#ifdef G4MULTITHREADED
  threadId = G4Threading::G4GetThreadId();
  #endif*/

  return m_map_writer[threadId];

}
  

bool GMCG4PodioManager::Finish()
{

  podio::ROOTWriter * l_writer = this->GetWriter();
  
  //  G4AutoLock lock(&GMCG4ActionMutex);

  if (l_writer != NULL) l_writer->finish();

  return true;
}

