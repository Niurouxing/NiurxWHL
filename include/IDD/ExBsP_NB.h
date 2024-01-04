#pragma once 
#include "IDDAlgorithm.h"
#include <vector>

class ExBsPCD;
class NBLDPC;
class MIMO;

class ExBsP_NB : public IDDAlgorithm<ExBsPCD, NBLDPC>
{
    using IDDAlgorithm<ExBsPCD, NBLDPC>::IDDAlgorithm;

};


