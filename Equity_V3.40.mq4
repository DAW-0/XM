//+------------------------------------------------------------------+
//|                                                       Equity.mq4 |
//|                          Copyright 2015-2025, JC LEGAZA Software |
//+------------------------------------------------------------------+

#property copyright   "2015-2025, JC LEGAZA Software"
#property description "Expert Advisor Equity v3.36"

#define Magic   2016011100

//+------------------------------------------------------------------+
//| Variables du programme                                           |
//+------------------------------------------------------------------+

extern int        TP    =12;        // Point de sortie
extern double     Risk  =0.9;       // Ratio de risque
extern bool       Trade =true;      // Ouvrir positions
extern bool       Hedge =true;      // Couvrir pertes
extern bool       Turbo =true;      // Profit anticipé
extern bool       Kill  =true;      // Solde balances

//  Tickets des positions ouvertes
//    Initialize the variables as arrays. MQL4 uses 0 indexing, but for readability we will leave index 0 as empty, i.e. BT0 = 0
int         BT[7] = {0, 0, 0, 0, 0, 0, 0}; // int        BT1, BT2, BT3, BT4, BT5, BT6;
int         ST[7] = {0, 0, 0, 0, 0, 0, 0}; // int        ST1, ST2, ST3, ST4, ST5, ST6;

//  Prix d'ouverture des positions
double      BP[7] = {0, 0, 0, 0, 0, 0, 0}; // double     BP1, BP2, BP3, BP4, BP5, BP6;
double      SP[7] = {0, 0, 0, 0, 0, 0, 0}; // double     SP1, SP2, SP3, SP4, SP5, SP6;

//  Gain des positions ouvertes
double      SG[7] = {0, 0, 0, 0, 0, 0, 0}; // double     SG1, SG2, SG3, SG4, SG5, SG6;
double      BG[7] = {0, 0, 0, 0, 0, 0, 0}; // double     BG1, BG2, BG3, BG4, BG5, BG6;

//  Volumes des positions ouvertes
double      BL[7] = {0, 0, 0, 0, 0, 0, 0}; // double     BL1, BL2, BL3, BL4, BL5, BL6;
double      SL[7] = {0, 0, 0, 0, 0, 0, 0}; // double     SL1, SL2, SL3, SL4, SL5, SL6;

//  Magic Number des positions ouvertes
int         BM[7] = {0, 0, 0, 0, 0, 0, 0}; // int        BM1, BM2, BM3, BM4, BM5, BM6;
int         SM[7] = {0, 0, 0, 0, 0, 0, 0}; // int        SM1, SM2, SM3, SM4, SM5, SM6;

//  Date d'ouverture des positions ouvertes
datetime    BO[7] = {0, 0, 0, 0, 0, 0, 0}; // datetime   BO1, BO2, BO3, BO4, BO5, BO6;
datetime    SO[7] = {0, 0, 0, 0, 0, 0, 0}; // datetime   SO1, SO2, SO3, SO4, SO5, SO6;

//  Paramètres système
double      Slip=30, Coef=1, Factor=1;
double      Cash=0, Equity=0;
double      BuyDD=0, SellDD=0, BuyDDtime=0, SellDDtime=0;
double      BuyTP=0, SellTP=0, OldBuyTP=0, OldSellTP=0;
double      HW, LW, CW, PP;
double      S[3] = {0, 0, 0};
double      R[3] = {0, 0, 0};
bool        Report=true, Success=false;
double      Pip=0, MinStep=0, Freeze=0;
datetime    LastTrade=0;
string      Status=" None";
int         Coin=1, Set=0;

//+------------------------------------------------------------------+
//| Check Status                                                     |
//+------------------------------------------------------------------+
void CheckStatus()
{
//    Reset Base de données Buy
//    2022-12-06 DAW ArrayInitialize is a built in function that will set any array to a chosen value. In this case, 0
      ArrayInitialize(BT, 0); // BT1=0; BT2=0; BT3=0; BT4=0; BT5=0; BT6=0;
      ArrayInitialize(BP, 0); // BP1=0; BP2=0; BP3=0; BP4=0; BP5=0; BP6=0;
      ArrayInitialize(BG, 0); // BG1=0; BG2=0; BG3=0; BG4=0; BG5=0; BG6=0;
      ArrayInitialize(BL, 0); // BL1=0; BL2=0; BL3=0; BL4=0; BL5=0; BL6=0;
      ArrayInitialize(BO, 0); // BO1=0; BO2=0; BO3=0; BO4=0; BO5=0; BO6=0;

//    Reset Base de données Sell
      ArrayInitialize(ST, 0); // ST1=0; ST2=0; ST3=0; ST4=0; ST5=0; ST6=0;
      ArrayInitialize(SP, 0); // SP1=0; SP2=0; SP3=0; SP4=0; SP5=0; SP6=0;
      ArrayInitialize(SG, 0); // SG1=0; SG2=0; SG3=0; SG4=0; SG5=0; SG6=0;
      ArrayInitialize(SL, 0); // SL1=0; SL2=0; SL3=0; SL4=0; SL5=0; SL6=0;
      ArrayInitialize(SO, 0); // SO1=0; SO2=0; SO3=0; SO4=0; SO5=0; SO6=0;
     
//    Recensement de tous les ordres ouverts
//    2022-12-06 DAW Now that the variables are arrays, we can loop through them instead of writing each one
      for(int i=OrdersTotal(); i>=1; i--) {
            if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==true) {
                  //  Renseigner la base de données des ordres Buy
                  BT[i] = OrderTicket(); BP[i] = OrderOpenPrice(); BL[i] = OrderLots(); BG[i] = MathFloor(OrderProfit() + OrderSwap()); BO[i] = OrderOpenTime();
                  //  Renseigner la base de données des ordres Sell
                  ST[i] = OrderTicket(); SP[i] = OrderOpenPrice(); SL[i] = OrderLots(); SG[i] = MathFloor(OrderProfit() + OrderSwap()); SO[i] = OrderOpenTime();
            }
      }
//    Calcul des DD Buy et Sell
      BuyDD = BG[1] + BG[2] + BG[3] + BG[4] + BG[5] + BG[6];
      SellDD= SG[1] + SG[2] + SG[3] + SG[4] + SG[5] + SG[6];
//    Calcul du Cash et de Equity
      Equity = MathRound(AccountEquity()-AccountCredit());
      Cash = MathRound(AccountBalance());
//    Calcul des High, Low et Close des bougies S et S-1 toutes les heures
      if(Minute()!=0 && Seconds()==0) {
            HW=MathMax(iHigh(Symbol(),PERIOD_W1,2),iHigh(Symbol(),PERIOD_W1,1));
            LW=MathMin(iLow(Symbol(),PERIOD_W1,2),iLow(Symbol(),PERIOD_W1,1));
            CW=(iClose(Symbol(),PERIOD_W1,2)+iClose(Symbol(),PERIOD_W1,1))/2;
//          Calcul du Point Pivot, Supports S1, S2, S3 et Résistances R1, R2, R3
            PP=(HW+LW+CW)/3;
            S[1]=2*PP-HW;
            S[2]=PP-(HW-LW);
            S[3]=PP-2*(HW-LW);
            R[1]=2*PP+LW;
            R[2]=PP+(HW-LW);
            R[3]=PP+2*(HW-LW);
//          Affichage des infos du compte et des positions
/*          if(IsTesting()==true)
                  Alert("Report Price ",Bid," Cash ",Cash," Equity ",Equity,
                  " BG1 ",BG1," BG2 ",BG2," BG3 ",BG3," BG4 ",BG4," BG5 ",BG5,
                  " BG6 ",BG6," SG1 ",SG1," SG2 ",SG2," SG3 ",SG3,
                  " SG4 ",SG4," SG5 ",SG5," SG6 ",SG6);
                  Report=false;
*/
      }
//    Taille du Pip
      Pip=10*Point;
//    Distance minimale autorisée pour placer un ordre
      MinStep=MarketInfo(Symbol(),MODE_STOPLEVEL)*Pip;
      Freeze=MarketInfo(Symbol(),MODE_FREEZELEVEL)*Pip;
//    Reporting de Cash et Equity
      Status=StringConcatenate(" Cash = ",Cash," Equity = ",Equity);
//    Si BuyTP-SellTP > 200 pips, Coef = BuyTP-SellTP / 200
      Coef=1;
      if(MathAbs(BuyTP-SellTP)>=0.02) {
            Coef=MathAbs(BuyTP-SellTP)/0.02;
      }
//    Calcul de Factor
      if(Coef!=0) {
            Factor=(Risk*Equity)/(1000*Coef);
      }
}
//+------------------------------------------------------------------+
//| Open Buy                                                         |
//+------------------------------------------------------------------+
void OpenBuy()
{
//  Si B1 n'existe pas
      CheckStatus();
      if(BP[1]==0) {
            BL[1]=NL(Factor*0.02);
            BT[1]=OrderSend(Symbol(),OP_BUY,BL[1],Ask,Slip,0,0,"B1",Magic+11,0,Green);
            if(IsTesting()==true) Alert("Open B1 = ",BL[1],Status);
      }
}
//+------------------------------------------------------------------+
//| Open Sell                                                        |
//+------------------------------------------------------------------+
void OpenSell()
{
//    Si S1 n'existe pas
      CheckStatus();
      if(SP[1]==0) {
            SL[1]=NL(Factor*0.02);
            ST[1]=OrderSend(Symbol(),OP_SELL,SL[1],Bid,Slip,0,0,"S1",Magic+21,0,Red);
            if(IsTesting()==true) Alert("Open S1 = ",SL[1],Status);
      }
}
//+------------------------------------------------------------------+
//| Trade Buy                                                        |
//+------------------------------------------------------------------+
void TradeBuy()
{
//        indices = {0,  1,  2,  3,  4,   5}
      int pips[6] = {0, 20, 30, 50, 80, 130};
      string txt = "";
//    Si B[i+1] n'existe pas et B[i] est en perte de pips[i], ouvrir B[i+1]
      CheckStatus();
      for(int i=OrdersTotal()-1; i>=1; i--) {
            if(BP[i+1]==0 && BP[i]!=0 && Ask<=BP[i]-Coef*pips[i]*10*Point) {
                  BL[i+1]=NL(Factor*pips[i]/1000);
                  if(IsTesting()==true) Alert("Trade B[",i+1,"] = ",BL[i+1],Status);
                  txt = StringConcatenate("B"+IntegerToString(i+1)); // make a string to pass to OrderSend
                  BT[i+1]=OrderSend(Symbol(),OP_BUY,BL[i+1],Ask,Slip,0,0,txt,Magic+11+i,0,Green);
            }
      }
}
//+------------------------------------------------------------------+
//| Trade Sell                                                       |
//+------------------------------------------------------------------+
void TradeSell()
{
//        indices = {0,  1,  2,  3,  4,   5}
      int pips[6] = {0, 20, 30, 50, 80, 130};
      string txt = "";
//    Si S[i+1] n'existe pas et S[i] est en perte de pips[i], ouvrir S[i+1]
      CheckStatus();
      for(int i=OrdersTotal()-1; i>=1; i--) {
            if(SP[i+1]==0 && SP[i]!=0 && Bid>=SP[i]+Coef*pips[i]*10*Point) {
                  SL[i+1]=NL(Factor*pips[i]/1000);
                  if(IsTesting()==true) Alert("Trade S[",i+1,"] = ",SL[i+1],Status);
                  txt = StringConcatenate("S"+IntegerToString(i+1)); // make a string to pass to OrderSend
                  ST[i+1]=OrderSend(Symbol(),OP_SELL,SL[i+1],Bid,Slip,0,0,txt,Magic+21+i,0,Red);
            }
      }

}
//+------------------------------------------------------------------+
//| Hedge Buy                                                        |
//+------------------------------------------------------------------+
void HedgeBuy()
{
//    indices = {0,  1,  2,  3,  4,   5}
      int pips[6] = {0, 30, 50, 70, 90, 110};
      string txt = "";
      if(Ask>=BP[1]) {
            for(int i=OrdersTotal()-1; i>=1; i--) {
      //          Si B[i+1] n'existe pas et B[i] est en gain de pips[i], ouvrir B2
                  CheckStatus();
                  if(BP[i]!=0 && BP[i+1]==0 && Ask>=BP[i]+pips[i]*10*Point) {
                        BL[i+1]=NL(BL[i]);
                        if(IsTesting()==true) Alert("Hedge B[",i+1,"] = ",BL[i+1],Status);
                        txt = StringConcatenate("B"+IntegerToString(i+1)); // make a string to pass to OrderSend
                        BT[i+1]=OrderSend(Symbol(),OP_BUY,BL[i+1],Ask,Slip,0,0,txt,Magic+11+i,0,Green);
                  }
            }
      }
}
//+------------------------------------------------------------------+
//| Hedge Sell                                                       |
//+------------------------------------------------------------------+
void HedgeSell()
{
//    indices = {0,  1,  2,  3,  4,   5}
      int pips[6] = {0, 30, 50, 70, 90, 110};
      string txt = "";
      if(Bid<=SP[1]) {
            for(int i=OrdersTotal()-1; i>=1; i--) {
//          Si S[i+1] n'existe pas et S[i] est en gain de pips[i], ouvrir S[i+1]
                  CheckStatus();
                  if(SP[i]!=0 && SP[i+1]==0 && Bid<=SP[i]-pips[i]*10*Point) {
                        SL[i+1]=NL(SL[i]);
                        if(IsTesting()==true) Alert("Hedge S[",i+1,"] = ",SL[i+1],Status);
                        txt = StringConcatenate("S"+IntegerToString(i+1)); // make a string to pass to OrderSend
                        ST[i+1]=OrderSend(Symbol(),OP_SELL,SL[i+1],Bid,Slip,0,0,txt,Magic+21+1,0,Red);
                  }
            }
      }
}
//+------------------------------------------------------------------+
//| Set TP Buy                                                       |
//+------------------------------------------------------------------+
void SetTPBuy()
{
//    2022-12-06 DAW I will need to make a summation and matrix multiplication function to streamline this. Skipping for now.
//    Calculer le TP commun à la moyenne des positions plus TP
      CheckStatus();
      if(BP[1]!=0 && BP[2]!=0) {
            BuyTP=NP(((BP[1]*BL[1]+BP[2]*BL[2]+BP[3]*BL[3]+BP[4]*BL[4]+BP[5]*BL[5]+BP[6]*BL[6])/(BL[1]+BL[2]+BL[3]+BL[4]+BL[5]+BL[6]))+TP*Pip);
      }
//    Modifier le TP commun des positions ouvertes
      if(OldBuyTP!=BuyTP &&
      MathAbs(BuyTP-Bid)>MinStep && MathAbs(BuyTP-Bid)>Freeze &&
      MathAbs(BP[1]-Bid)>MinStep && MathAbs(BP[1]-Bid)>Freeze &&
      MathAbs(BP[2]-Bid)>MinStep && MathAbs(BP[2]-Bid)>Freeze &&
      MathAbs(BP[3]-Bid)>MinStep && MathAbs(BP[3]-Bid)>Freeze &&
      MathAbs(BP[4]-Bid)>MinStep && MathAbs(BP[4]-Bid)>Freeze &&
      MathAbs(BP[5]-Bid)>MinStep && MathAbs(BP[5]-Bid)>Freeze &&
      MathAbs(BP[6]-Bid)>MinStep && MathAbs(BP[6]-Bid)>Freeze) {
            for(int i=OrdersTotal(); i>=1; i--) {
                  if(BP[i]!=0 && Bid<BuyTP) Success=OrderModify(BT[i],0,0,BuyTP,0,Green);
                  if(BP[i]!=0 && Bid>BuyTP) Success=OrderModify(BT[i],0,BuyTP,0,0,Green);
            }
//          Enregistrer le TP commun des positions ouvertes
            if(IsTesting()==true) Alert("Set BuyTP = ",BuyTP);
            LastTrade=TimeCurrent();
            OldBuyTP=BuyTP;
      }
}
//+------------------------------------------------------------------+
//| Set TP Sell                                                      |
//+------------------------------------------------------------------+
void SetTPSell()
{
//    2022-12-06 DAW I will need to make a summation and matrix multiplication function to streamline this. Skipping for now.
//    Calculer le TP commun à la moyenne des positions moins TP
      CheckStatus();
      if(SP[1]!=0 && SP[2]!=0) {
            SellTP=NP(((SP[1]*SL[1]+SP[2]*SL[2]+SP[3]*SL[3]+SP[4]*SL[4]+SP[5]*SL[5]+SP[6]*SL[6])/(SL[1]+SL[2]+SL[3]+SL[4]+SL[5]+SL[6]))-TP*Pip);
      }
//    Modifier le TP commun des positions ouvertes
      if(OldSellTP!=SellTP &&
      MathAbs(SellTP-Ask)>MinStep &&
      MathAbs(SellTP-Ask)>Freeze  &&
      MathAbs(SP[1] -Ask)>MinStep &&
      MathAbs(SP[2] -Ask)>MinStep &&
      MathAbs(SP[3] -Ask)>MinStep &&
      MathAbs(SP[4] -Ask)>MinStep &&
      MathAbs(SP[5] -Ask)>MinStep &&
      MathAbs(SP[6] -Ask)>MinStep) {
            for(int i=OrdersTotal(); i>=1; i--) {
                  if(SP[i]!=0 && Ask>SellTP) Success=OrderModify(ST[i],0,0,SellTP,0,Red);
                  if(SP[i]!=0 && Ask<SellTP) Success=OrderModify(ST[i],0,SellTP,0,0,Red);
            }
//          Enregistrer le TP commun des positions ouvertes
            if(IsTesting()==true) Alert("Set SellTP = ",SellTP);
            LastTrade=TimeCurrent();
            OldSellTP=SellTP;
      }
}
//+------------------------------------------------------------------+
//| Turbo Buy                                                        |
//+------------------------------------------------------------------+
void TurboBuy()
{
//    2022-12-06 DAW I will need to make a summation function to streamline this. Skipping for now.
//    Si BuyDD < 35% Risk*Cash
      if(BuyDD>-0.35*Risk*Cash) {
//    Si B3 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(BP[3]!=0 && BG[3]>=(BL[1]+BL[2]+BL[3])*10*TP && TimeCurrent()-BO[2]<=3*24*3600) {
            if(IsTesting()==true) Alert("Turbo B3 = ",BG[3]);
            Success=OrderClose(BT[3],BL[3],Bid,Slip,Green);
      }
//    Si B4 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(BP[4]!=0 && BG[4]>=(BL[1]+BL[2]+BL[3]+BL[4])*10*TP && TimeCurrent()-BO[3]<=5*24*3600) {
            if(IsTesting()==true) Alert("Turbo B4 = ",BG[4]);
            Success=OrderClose(BT[4],BL[4],Bid,Slip,Green);
      }
//    Si B5 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(BP[5]!=0 && BG[5]>=(BL[1]+BL[2]+BL[3]+BL[4]+BL[5])*10*TP && TimeCurrent()-BO[4]<=8*24*3600) {
      if(IsTesting()==true) Alert("Turbo B5 = ",BG[5]);
      Success=OrderClose(BT[5],BL[5],Bid,Slip,Green);
      }
      }
   }
//+------------------------------------------------------------------+
//| Turbo Sell                                                       |
//+------------------------------------------------------------------+
void TurboSell()
{
//    2022-12-06 DAW I will need to make a summation function to streamline this. Skipping for now.
//    Si SellDD < 35% Risk*Cash
      if(SellDD>-0.35*Risk*Cash) {
//          Si S3 fait le gain du TP commun, la cloturer
            CheckStatus();
            if(SP[3]!=0 && SG[3]>=(SL[1]+SL[2]+SL[3])*10*TP && TimeCurrent()-SO[2]<=3*24*3600) {
                  Success=OrderClose(ST[3],SL[3],Ask,Slip,Red);
                  if(IsTesting()==true) Alert("Turbo S3 = ",SG[3]);
            }
//          Si S4 fait le gain du TP commun, la cloturer
            CheckStatus();
            if(SP[4]!=0 && SG[4]>=(SL[1]+SL[2]+SL[3]+SL[4])*10*TP && TimeCurrent()-SO[3]<=5*24*3600) {
                  Success=OrderClose(ST[4],SL[4],Ask,Slip,Red);
                  if(IsTesting()==true) Alert("Turbo S4 = ",SG[4]);
            }
//          Si S5 fait le gain du TP commun, la cloturer
            CheckStatus();
            if(SP[5]!=0 && SG[5]>=(SL[1]+SL[2]+SL[3]+SL[4]+SL[5])*10*TP && TimeCurrent()-SO[4]<=8*24*3600) {
                  Success=OrderClose(ST[5],SL[5],Ask,Slip,Red);
                  if(IsTesting()==true) Alert("Turbo S5 = ",SG[5]);
            }
      }
}
//+------------------------------------------------------------------+
//| Close TP Buy                                                     |
//+------------------------------------------------------------------+
void CloseTPBuy()
{
//    2022-12-06 DAW Really need to write a matrix multiplication function as well to clean this up... Grrrrrrr.... MQL is so basic :(
//    Calculer le TP commun à la moyenne des positions plus TP
      CheckStatus();
      if(BL[1]!=0){
            BuyTP=NP(((BP[1]*BL[1]+BP[2]*BL[2]+BP[3]*BL[3]+BP[4]*BL[4]+BP[5]*BL[5]+BP[6]*BL[6])/(BL[1]+BL[2]+BL[3]+BL[4]+BL[5]+BL[6]))+TP*10*Point);
      }
//    Clôturer les positions au TP commun
      if(BL[2]!=0 && Bid>=BuyTP) {
            if(IsTesting()==true) Alert("Close TP Buy = ",BuyDD,Status);
            CloseBuy();
      }
}
//+------------------------------------------------------------------+
//| Close TP Sell                                                    |
//+------------------------------------------------------------------+
void CloseTPSell()
{
//    Calculer le TP commun à la moyenne des positions moins TP
      CheckStatus();
      if(SL[1]!=0){
            SellTP=NP(((SP[1]*SL[1]+SP[2]*SL[2]+SP[3]*SL[3]+SP[4]*SL[4]+SP[5]*SL[5]+SP[6]*SL[6])/(SL[1]+SL[2]+SL[3]+SL[4]+SL[5]+SL[6]))-TP*10*Point);
      }
//    Clôturer les positions au TP commun
      if(SL[2]!=0 && Ask<=SellTP) {
            if(IsTesting()==true) Alert("Close TP Sell = ",SellDD,Status);
            CloseSell();
      }
}
//+------------------------------------------------------------------+
//| Close Buy                                                        |
//+------------------------------------------------------------------+
void CloseBuy()
{
//    Clôture des Buy
      CheckStatus();
      for(int i=OrdersTotal(); i>=1; i--) {
            if(BP[i]!=0) Success=OrderClose(BT[i],BL[i],Bid,Slip,Green);
      }
      if(IsTesting()==true) Alert("CloseBuy = ",BG[1]+BG[2]+BG[3]+BG[4]+BG[5]+BG[6]);
}
//+------------------------------------------------------------------+
//| Close Sell                                                       |
//+------------------------------------------------------------------+
void CloseSell()
{
//    Clôture des Sell
      CheckStatus();
      for(int i=OrdersTotal(); i>=1; i--) {
            if(SP[i]!=0) Success=OrderClose(ST[i],SL[i],Ask,Slip,Red);
      }
      if(IsTesting()==true) Alert("CloseSell = ",SG[1]+SG[2]+SG[3]+SG[4]+SG[5]+SG[6]);
}
//+------------------------------------------------------------------+
//| Close Buy                                                        |
//+------------------------------------------------------------------+
void CloseBuyOld()
{
//    Calculer le TP commun à la moyenne des positions ponérées plus TP
      CheckStatus();
      if(BP[1]!=0) {
            BuyTP=NP(((BP[1]*BL[1]+BP[2]*BL[2]+BP[3]*BL[3]+BP[4]*BL[4]+BP[5]*BL[5]+BP[6]*BL[6])/(BL[1]+BL[2]+BL[3]+BL[4]+BL[5]+BL[6]))+TP*10*Point);
      }
//    Modifier le TP commun des positions ouvertes
      if(Bid>=BuyTP) {
            for(int i=OrdersTotal(); i>=1; i--) {
                  if(BP[i]!=0) Success=OrderClose(BT[i],BL[i],Bid,Slip,Green);
            }
            if(IsTesting()==true) Alert("CloseBuy = ",BG[1]+BG[2]+BG[3]+BG[4]+BG[5]+BG[6]);
      }
}
//+------------------------------------------------------------------+
//| Close Sell                                                       |
//+------------------------------------------------------------------+
void CloseSellOld()
{
//    Calculer le TP commun à ma moyenne des positions ponérées moins TP
      CheckStatus();
      if(SP[1]!=0) {
            SellTP=NP(((SP[1]*SL[1]+SP[2]*SL[2]+SP[3]*SL[3]+SP[4]*SL[4]+SP[5]*SL[5]+SP[6]*SL[6])/(SL[1]+SL[2]+SL[3]+SL[4]+SL[5]+SL[6]))-TP*10*Point);
      }
//    Modifier le TP commun des positions ouvertes
      if(Ask<=SellTP) {
            for(int i=OrdersTotal(); i>=1; i--) {
                  if(SP[i]!=0) Success=OrderClose(ST[i],SL[i],Ask,Slip,Red);
            }
            if(IsTesting()==true) Alert("CloseSell = ",SG[1]+SG[2]+SG[3]+SG[4]+SG[5]+SG[6]);
      }
}
//+------------------------------------------------------------------+
//| Normalize Price                                                  |
//+------------------------------------------------------------------+
double NP(double price)
{
//    Normalisation du prix à la cinquième décimale
      price=MathFloor(price*100000)/100000;
      return (price);
}
//+------------------------------------------------------------------+
//| Normalize Lot                                                    |
//+------------------------------------------------------------------+
double NL(double lot)
{   
//    Normalisation du lot à la deuxième décimale, lot mini 0.01
      lot=MathFloor(lot*100)/100;
      if(lot<=0.02) lot=0.02;
      return (lot);
}
//+------------------------------------------------------------------+
//| Information                                                      |
//+------------------------------------------------------------------+
void Info()
{
      string BarName;
      
//    Reset des infos sur le graphe
      ObjectsDeleteAll();
//    Position des infos sur le graphe
      if(Ask>(WindowPriceMax()+WindowPriceMin())/1.9) Coin=1;
      if(Ask<(WindowPriceMax()+WindowPriceMin())/2.1) Coin=3;

      CheckStatus();
      BarName = "Parameters1";
      ObjectCreate(BarName,OBJ_LABEL,0,0,0);
      ObjectSetText(BarName,"Risk="+NL(Risk)+" TP="+TP,12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,400);
      ObjectSet(BarName,OBJPROP_YDISTANCE,20);

      CheckStatus();
      BarName = "Parameters2";
      ObjectCreate(BarName,OBJ_LABEL,0,0,0);
      ObjectSetText(BarName,"Coef="+NL(Coef)+" Factor="+NL(Factor),12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,400);
      ObjectSet(BarName,OBJPROP_YDISTANCE,40);

      BarName = "Volumes";
      ObjectCreate(BarName,OBJ_LABEL,0,0,0);
      ObjectSetText(BarName,"BuyTP= "+BuyTP+" SellTP= "+SellTP,12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,400);
      ObjectSet(BarName,OBJPROP_YDISTANCE,60);

      BarName = "DD";
      ObjectCreate(BarName,OBJ_LABEL,0,0,0);
      ObjectSetText(BarName,"BDD= "+MathRound(100*(BuyDD)/Cash)+" %"+" SDD = "+MathRound(100*(SellDD)/Cash)+" %",12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,400);
      ObjectSet(BarName,OBJPROP_YDISTANCE,80);

      CheckStatus();
      BarName = "Money";
      ObjectCreate(BarName,OBJ_LABEL,0,0,0);
      ObjectSetText(BarName,"Cash= "+MathRound(Cash)+" Equity= "+MathRound(AccountEquity()),12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,400);
      ObjectSet(BarName,OBJPROP_YDISTANCE,100);

      BarName = "PP";
      ObjectDelete(BarName);
      ObjectCreate(BarName,OBJ_HLINE,0,0,PP);
      ObjectSet(BarName,OBJPROP_COLOR,Blue);
      ObjectSet(BarName,OBJPROP_WIDTH,1);
      ObjectSet(BarName,OBJPROP_STYLE,3);

      string aTxt[3] = {"", "R", "S"};
      string txt = "";
      for(int i=3; i>=1; i--) {
            for(int j=2; j>=1; j--) {
                  txt = StringConcatenate(aTxt[j]+IntegerToString(i));
                  BarName = txt;
                  ObjectDelete(BarName);
                  if(j==1) {
                        ObjectCreate(BarName,OBJ_HLINE,0,0,R[i]);
                  }
                  else if(j==2) {
                        ObjectCreate(BarName,OBJ_HLINE,0,0,S[i]);
                  }
                  ObjectSet(BarName,OBJPROP_COLOR,Yellow);
                  ObjectSet(BarName,OBJPROP_WIDTH,1);
                  ObjectSet(BarName,OBJPROP_STYLE,3);
            }
      }
}
//+------------------------------------------------------------------+
//| Execute                                                          |
//+------------------------------------------------------------------+
void OnTick()
{
//    Trader à Minute=0 et Seconde=0 
      if((IsTesting()==true && Seconds()==0) || (IsTesting()==false && IsTradeAllowed()==true)) {
            OpenBuy();
            OpenSell();
//          Trader à Minute=0 et Seconde=0 et si Equity > 50
            if(Trade==true && Minute()==0 && Seconds()==0 && Equity>=50)
//          Si le marché est calme (volatilité M5 < TP pips) depuis 5 minutes
            if(iHigh(Symbol(),PERIOD_M5,1)-iLow(Symbol(),PERIOD_M5,1)<=(TP*10*Point))
//          S'il est moins de 20h
            if(Hour()<20 && Hour()>8) {
                  TradeBuy();
                  TradeSell();
            }
//          Mettre à jour le TP de sortie des martingales
            SetTPBuy();
            SetTPSell();
            CloseTPBuy();
            CloseTPSell();
//          Couvrir la martingale perdante par une gagnante 
            if(Hedge==true) {
            HedgeBuy();
            HedgeSell();
            }
//          Prendre le gain de la position gagnante
            if(Turbo==true) {
            TurboBuy();
            TurboSell();
            }
//          Clôturer les DD si gain=pertes
            if(Kill==true && Equity>Cash && (BuyDD<-0.45*Risk*Cash || SellDD<-0.45*Risk*Cash)) {
            CloseBuy();
            CloseSell();
            }
//          Afficher les infos des martingales
            Info();
      }
}