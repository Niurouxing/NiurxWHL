#include "ExBsP_NB.h"
#include "ExBsP.h"
#include "NBLDPC.h"
#include "MIMO.h"


// void ExBsP_NB::bind(MIMO * mimo)
// {
//     nbldpc = new NBLDPC();
//     nbldpc->bind(mimo);
//     for(int i=0;i<mimo->getN();i++){
//         ExBsPCD *exBspcd = new ExBsPCD();
//         exBspcd->bind(mimo->detections[i]);
//         exBspcds.push_back(exBspcd);
//     }
// }