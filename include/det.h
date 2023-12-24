#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Detection.h"
#include "ExBsP.h"


Detection* detection;
DetectionAlgorithm* algorithm;


void detectionInit(int TxAntNum, int RxAntNum, int ModType, double SNRdB, bool isComplex){
    if(isComplex){
        detection = new DetectionCD(TxAntNum, RxAntNum, ModType, SNRdB);
    }else{
        detection = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
    }
}

void MMSEInit(bool isComplex){
    if(isComplex){
        algorithm = new MMSECD();
    }else{
        algorithm = new MMSE();
    }
    algorithm->bind(detection);
}

void EPInit(int iter, double delta){
    algorithm = new EP(iter, delta);
    algorithm->bind(detection);
}

void ExBsPInit(int iter, int dm){
    algorithm = new ExBsPCD(iter, dm);
    algorithm->bind(detection);
}

void execute(){
    detection->generate();
    algorithm->execute();
    algorithm->check();
}

std::tuple<int,int> report(){
    int newErrorBits = algorithm->getNewErrorBits();
    int newErrorFrames = algorithm->getNewErrorFrames();

    return std::make_tuple(newErrorBits, newErrorFrames);
}

 