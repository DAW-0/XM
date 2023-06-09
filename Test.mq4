//+------------------------------------------------------------------+
//|                                                        Test1.mq4 |
//|                        Copyright 2022, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property strict

//  Tickets des positions ouvertes
int      BT[6] = {0, 0, 0, 0, 0, 0};

//+------------------------------------------------------------------+
//| Expert initialization function                                   |
//+------------------------------------------------------------------+
int OnInit()
  {
//---
//   for( i=3; i>=0; i--){
    Print("Testing Init 1 2 3");
    Print("Initialized BT =", BT[1]);
    BT[1] = 999;
    Print("Set BT[1] =", BT[1]);
    ArrayInitialize(BT, 0); // BT1=0; BT2=0; BT3=0; BT4=0; BT5=0; BT6=0;
    Print("Reinitialized BT =", BT[1]);
//   }
//---
   return(INIT_SUCCEEDED);
  }
//+------------------------------------------------------------------+
//| Expert deinitialization function                                 |
//+------------------------------------------------------------------+
void OnDeinit(const int reason)
  {
//---
   Print("Testing Deinit 1 2 3");
  }
//+------------------------------------------------------------------+
//| Expert tick function                                             |
//+------------------------------------------------------------------+
void OnTick()
  {
//---
   Print("Testing Ontick 1 2 3");
  }
//+------------------------------------------------------------------+