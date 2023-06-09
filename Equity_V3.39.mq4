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

extern int      TP          =12;          // Point de sortie
extern double   Risk        =0.9;         // Ratio de risque
extern bool     Trade       =true;        // Ouvrir positions
extern bool     Hedge       =true;        // Couvrir pertes
extern bool     Turbo       =true;        // Profit anticipé
extern bool     Kill        =true;        // Solde balances

//  Tickets des positions ouvertes
int      BT1, BT2, BT3, BT4, BT5, BT6;
int      ST1, ST2, ST3, ST4, ST5, ST6;

//  Prix d'ouverture des positions
double   BP1, BP2, BP3, BP4, BP5, BP6;
double   SP1, SP2, SP3, SP4, SP5, SP6;

//  Gain des positions ouvertes
double   BG1, BG2, BG3, BG4, BG5, BG6;
double   SG1, SG2, SG3, SG4, SG5, SG6;

//  Volumes des positions ouvertes
double   BL1, BL2, BL3, BL4, BL5, BL6;
double   SL1, SL2, SL3, SL4, SL5, SL6;

//  Magic Number des positions ouvertes
int      BM1, BM2, BM3, BM4, BM5, BM6;
int      SM1, SM2, SM3, SM4, SM5, SM6;

//  Date d'ouverture des positions ouvertes
datetime BO1, BO2, BO3, BO4, BO5, BO6;
datetime SO1, SO2, SO3, SO4, SO5, SO6;

//  Paramètres système
double   Slip=30, Coef=1, Factor=1;
double   Cash=0, Equity=0;
double   BuyDD=0, SellDD=0, BuyDDtime=0, SellDDtime=0;
double   BuyTP=0, SellTP=0, OldBuyTP=0, OldSellTP=0;
double   HW, LW, CW, PP, S1, S2, S3, R1, R2, R3;
bool     Report=true, Success=false;
double   Pip=0, MinStep=0, Freeze=0;
datetime LastTrade=0;
string   Status=" None";
int      Coin=1, Set=0;

//+------------------------------------------------------------------+
//| Check Status                                                     |
//+------------------------------------------------------------------+
void CheckStatus()
  {
//  Reset Base de données Buy
      BT1=0; BT2=0; BT3=0; BT4=0; BT5=0; BT6=0;
      BP1=0; BP2=0; BP3=0; BP4=0; BP5=0; BP6=0;
      BG1=0; BG2=0; BG3=0; BG4=0; BG5=0; BG6=0;
      BL1=0; BL2=0; BL3=0; BL4=0; BL5=0; BL6=0;
      BO1=0; BO2=0; BO3=0; BO4=0; BO5=0; BO6=0;

//  Reset Base de données Sell
      ST1=0; ST2=0; ST3=0; ST4=0; ST5=0; ST6=0;
      SP1=0; SP2=0; SP3=0; SP4=0; SP5=0; SP6=0;
      SG1=0; SG2=0; SG3=0; SG4=0; SG5=0; SG6=0;
      SL1=0; SL2=0; SL3=0; SL4=0; SL5=0; SL6=0;
      SO1=0; SO2=0; SO3=0; SO4=0; SO5=0; SO6=0;
     
//  Recensement de tous les ordres ouverts
   for(int i=OrdersTotal()-1; i>=0; i--)
      if(OrderSelect(i,SELECT_BY_POS,MODE_TRADES)==true) {
//  Renseigner la base de données des ordres Buy
      if(OrderComment()=="B1" || OrderMagicNumber()==Magic+11) {
      BT1=OrderTicket(); BP1=OrderOpenPrice(); BL1=OrderLots(); BG1=MathFloor(OrderProfit()+OrderSwap()); BO1=OrderOpenTime(); }
      if(OrderComment()=="B2" || OrderMagicNumber()==Magic+12) {
      BT2=OrderTicket(); BP2=OrderOpenPrice(); BL2=OrderLots(); BG2=MathFloor(OrderProfit()+OrderSwap()); BO2=OrderOpenTime(); }
      if(OrderComment()=="B3" || OrderMagicNumber()==Magic+13) {
      BT3=OrderTicket(); BP3=OrderOpenPrice(); BL3=OrderLots(); BG3=MathFloor(OrderProfit()+OrderSwap()); BO3=OrderOpenTime(); }
      if(OrderComment()=="B4" || OrderMagicNumber()==Magic+14) {
      BT4=OrderTicket(); BP4=OrderOpenPrice(); BL4=OrderLots(); BG4=MathFloor(OrderProfit()+OrderSwap()); BO4=OrderOpenTime(); }
      if(OrderComment()=="B5" || OrderMagicNumber()==Magic+15) {
      BT5=OrderTicket(); BP5=OrderOpenPrice(); BL5=OrderLots(); BG5=MathFloor(OrderProfit()+OrderSwap()); BO5=OrderOpenTime(); }
      if(OrderComment()=="B6" || OrderMagicNumber()==Magic+16) {
      BT6=OrderTicket(); BP6=OrderOpenPrice(); BL6=OrderLots(); BG6=MathFloor(OrderProfit()+OrderSwap()); BO6=OrderOpenTime(); }
//  Renseigner la base de données des ordres Sell
      if(OrderComment()=="S1" || OrderMagicNumber()==Magic+21) {
      ST1=OrderTicket(); SP1=OrderOpenPrice(); SL1=OrderLots(); SG1=MathFloor(OrderProfit()+OrderSwap()); SO1=OrderOpenTime(); }
      if(OrderComment()=="S2" || OrderMagicNumber()==Magic+22) {
      ST2=OrderTicket(); SP2=OrderOpenPrice(); SL2=OrderLots(); SG2=MathFloor(OrderProfit()+OrderSwap()); SO2=OrderOpenTime(); }
      if(OrderComment()=="S3" || OrderMagicNumber()==Magic+23) {
      ST3=OrderTicket(); SP3=OrderOpenPrice(); SL3=OrderLots(); SG3=MathFloor(OrderProfit()+OrderSwap()); SO3=OrderOpenTime(); }
      if(OrderComment()=="S4" || OrderMagicNumber()==Magic+24) {
      ST4=OrderTicket(); SP4=OrderOpenPrice(); SL4=OrderLots(); SG4=MathFloor(OrderProfit()+OrderSwap()); SO4=OrderOpenTime(); }
      if(OrderComment()=="S5" || OrderMagicNumber()==Magic+25) {
      ST5=OrderTicket(); SP5=OrderOpenPrice(); SL5=OrderLots(); SG5=MathFloor(OrderProfit()+OrderSwap()); SO5=OrderOpenTime(); }
      if(OrderComment()=="S6" || OrderMagicNumber()==Magic+26) {
      ST6=OrderTicket(); SP6=OrderOpenPrice(); SL6=OrderLots(); SG6=MathFloor(OrderProfit()+OrderSwap()); SO6=OrderOpenTime(); }
      }
//  Calcul des DD Buy et Sell
      BuyDD=BG1+BG2+BG3+BG4+BG5+BG6;
      SellDD=SG1+SG2+SG3+SG4+SG5+SG6;
//  Calcul du Cash et de Equity
      Equity=MathRound(AccountEquity()-AccountCredit());
      Cash=MathRound(AccountBalance());
//  Calcul des High, Low et Close des bougies S et S-1 toutes les heures
      if(Minute()!=0 && Seconds()==0) Report=true;
      if(Report==true) {
      HW=MathMax(iHigh(Symbol(),PERIOD_W1,2),iHigh(Symbol(),PERIOD_W1,1));
      LW=MathMin(iLow(Symbol(),PERIOD_W1,2),iLow(Symbol(),PERIOD_W1,1));
      CW=(iClose(Symbol(),PERIOD_W1,2)+iClose(Symbol(),PERIOD_W1,1))/2;
//  Calcul du Point Pivot, Supports S1, S2, S3 et Résistances R1, R2, R3
      PP=(HW+LW+CW)/3;
      S1=2*PP-HW;
      S2=PP-(HW-LW);
      S3=PP-2*(HW-LW);
      R1=2*PP+LW;
      R2=PP+(HW-LW);
      R3=PP+2*(HW-LW);
//  Affichage des infos du compte et des positions
/*      if(IsTesting()==true)
      Alert("Report Price ",Bid," Cash ",Cash," Equity ",Equity,
      " BG1 ",BG1," BG2 ",BG2," BG3 ",BG3," BG4 ",BG4," BG5 ",BG5,
      " BG6 ",BG6," SG1 ",SG1," SG2 ",SG2," SG3 ",SG3,
      " SG4 ",SG4," SG5 ",SG5," SG6 ",SG6);
      Report=false;
*/
      }
//  Taille du Pip
      Pip=10*Point;
//  Distance minimale autorisée pour placer un ordre
      MinStep=MarketInfo(Symbol(),MODE_STOPLEVEL)*Pip;
      Freeze=MarketInfo(Symbol(),MODE_FREEZELEVEL)*Pip;
//  Reporting de Cash et Equity
      Status=StringConcatenate(" Cash = ",Cash," Equity = ",Equity);
//  Si BuyTP-SellTP > 200 pips, Coef = BuyTP-SellTP / 200
      Coef=1;
      if(MathAbs(BuyTP-SellTP)>=0.02)
      Coef=MathAbs(BuyTP-SellTP)/0.02;
//  Calcul de Factor
      if(Coef!=0)
      Factor=(Risk*Equity)/(1000*Coef);
  }
//+------------------------------------------------------------------+
//| Open Buy                                                         |
//+------------------------------------------------------------------+
void OpenBuy()
   {
//  Si B1 n'existe pas
      CheckStatus();
      if(BP1==0) {
      BL1=NL(Factor*0.02);
      BT1=OrderSend(Symbol(),OP_BUY,BL1,Ask,Slip,0,0,"B1",Magic+11,0,Green);
      if(IsTesting()==true) Alert("Open B1 = ",BL1,Status);
      }
   }
//+------------------------------------------------------------------+
//| Open Sell                                                        |
//+------------------------------------------------------------------+
void OpenSell()
   {
//  Si S1 n'existe pas
      CheckStatus();
      if(SP1==0) {
      SL1=NL(Factor*0.02);
      ST1=OrderSend(Symbol(),OP_SELL,SL1,Bid,Slip,0,0,"S1",Magic+21,0,Red);
      if(IsTesting()==true) Alert("Open S1 = ",SL1,Status);
      }
   }
//+------------------------------------------------------------------+
//| Trade Buy                                                        |
//+------------------------------------------------------------------+
void TradeBuy()
   { 
//  Si B2 n'existe pas et B1 est en perte de 20 pips, ouvrir B2
      CheckStatus();
      if(BP2==0 && BP1!=0 && Ask<=BP1-Coef*20*10*Point) {
      BL2=NL(Factor*0.02);
      if(IsTesting()==true) Alert("Trade B2 = ",BL2,Status);
      BT2=OrderSend(Symbol(),OP_BUY,BL2,Ask,Slip,0,0,"B2",Magic+12,0,Green);
      }
//  Si B3 n'existe pas et B2 est en perte de 30 pips, ouvrir B3
      CheckStatus();
      if(BP3==0 && BP2!=0 && Ask<=BP2-Coef*30*10*Point) {
      BL3=NL(Factor*0.03);
      if(IsTesting()==true) Alert("Trade B3 = ",BL3,Status);
      BT3=OrderSend(Symbol(),OP_BUY,BL3,Ask,Slip,0,0,"B3",Magic+13,0,Green);
      }
//  Si B4 n'existe pas, que S3 existe et B3 est en perte de 50 pips, ouvrir B4
      CheckStatus();
      if(BP4==0 && SP3!=0 && BP3!=0 && Ask<=BP3-Coef*50*10*Point) {
      BL4=NL(Factor*0.05);
      if(IsTesting()==true) Alert("Trade B4 = ",BL4,Status);
      BT4=OrderSend(Symbol(),OP_BUY,BL4,Ask,Slip,0,0,"B4",Magic+14,0,Green);
      }
//  Si B5 n'existe pas et B4 est en perte de 80 pips, ouvrir B5
      CheckStatus();
      if(BP5==0 && BP4!=0 && Ask<=BP4-Coef*80*10*Point) {
      BL5=NL(Factor*0.08);
      if(IsTesting()==true) Alert("Trade B5 = ",BL5,Status);
      BT5=OrderSend(Symbol(),OP_BUY,BL5,Ask,Slip,0,0,"B5",Magic+15,0,Green);
      }
//  Si B6 n'existe pas et B5 est en perte de 130 pips, ouvrir B6
      CheckStatus();
      if(BP6==0 && BP5!=0 && Ask<=BP5-Coef*130*10*Point) {
      if(IsTesting()==true) Alert("Trade B6 = ",BL6,Status);
      BL6=NL(Factor*0.13);
      BT6=OrderSend(Symbol(),OP_BUY,BL6,Ask,Slip,0,0,"B6",Magic+16,0,Green);
      }
   }
//+------------------------------------------------------------------+
//| Trade Sell                                                       |
//+------------------------------------------------------------------+
void TradeSell()
   {
//  Si S2 n'existe pas et S1 est en perte de 20 pips, ouvrir S2
      CheckStatus();
      if(SP2==0 && SP1!=0 && Bid>=SP1+Coef*20*10*Point) {
      SL2=NL(Factor*0.02);
      if(IsTesting()==true) Alert("Trade S2 = ",SL2,Status);
      ST2=OrderSend(Symbol(),OP_SELL,SL2,Bid,Slip,0,0,"S2",Magic+22,0,Red);
      }
//  Si S3 n'existe pas et S2 est en perte de 30 pips, ouvrir S3
      CheckStatus();
      if(SP3==0 && SP2!=0 && Bid>=SP2+Coef*30*10*Point) {
      SL3=NL(Factor*0.03);
      if(IsTesting()==true) Alert("Trade S3 = ",SL3,Status);
      ST3=OrderSend(Symbol(),OP_SELL,SL3,Bid,Slip,0,0,"S3",Magic+23,0,Red);
      }
//  Si S4 n'existe pas, que S3 existe et B3 est en perte de 50 pips, ouvrir S4
      CheckStatus();
      if(SP4==0 && SP3!=0 && BP3!=0 && Bid>=SP3+Coef*50*10*Point) {
      SL4=NL(Factor*0.05);
      if(IsTesting()==true) Alert("Trade S4 = ",SL4,Status);
      ST4=OrderSend(Symbol(),OP_SELL,SL4,Bid,Slip,0,0,"S4",Magic+24,0,Red);
      }
//  Si S5 n'existe pas et S4 est en perte de 80 pips, ouvrir S5
      CheckStatus();
      if(SP5==0 && SP4!=0 && Bid>=SP4+Coef*80*10*Point) {
      SL5=NL(Factor*0.08);
      if(IsTesting()==true) Alert("Trade S5 = ",SL5,Status);
      ST5=OrderSend(Symbol(),OP_SELL,SL5,Bid,Slip,0,0,"S5",Magic+25,0,Red);
      }
//  Si S6 n'existe pas et S5 est en perte de 130 pips, ouvrir S6
      CheckStatus();
      if(SP6==0 && SP5!=0 && Bid>=SP5+Coef*130*10*Point) {
      SL6=NL(Factor*0.13);
      if(IsTesting()==true) Alert("Trade S6 = ",SL6,Status);
      ST6=OrderSend(Symbol(),OP_SELL,SL6,Bid,Slip,0,0,"S6",Magic+26,0,Red);
      }
   }
//+------------------------------------------------------------------+
//| Hedge Buy                                                        |
//+------------------------------------------------------------------+
void HedgeBuy()
   { 
      if(Ask>=BP1) {
//  Si B2 n'existe pas et B1 est en gain de 30 pips, ouvrir B2
      CheckStatus();
      if(BP1!=0 && BP2==0 && Ask>=BP1+30*10*Point) {
      BL2=NL(BL1);
      if(IsTesting()==true) Alert("Hedge B2 = ",BL2,Status);
      BT2=OrderSend(Symbol(),OP_BUY,BL2,Ask,Slip,0,0,"B2",Magic+12,0,Green);
      }
//  Si B3 n'existe pas et B1 est en gain de 50 pips, ouvrir B3
      if(BP2!=0 && BP3==0 && Ask>=BP1+50*10*Point) {
      BL3=NL(SL1+SL2+SL3+SL4+SL5+SL6-BL2);
      if(IsTesting()==true) Alert("Hedge B3 = ",BL3,Status);
      BT3=OrderSend(Symbol(),OP_BUY,BL3,Ask,Slip,0,0,"B3",Magic+13,0,Green);
      }
//  Si B4 n'existe pas et B1 est en gain de 70 pips, ouvrir B4
      CheckStatus();
      if(BP3!=0 && BP4==0 && Ask>=BP1+70*10*Point) {
      BL4=NL(SL1+SL2+SL3+SL4+SL5+SL6-BL3);
      if(IsTesting()==true) Alert("Hedge B4 = ",BL4,Status);
      BT4=OrderSend(Symbol(),OP_BUY,BL4,Ask,Slip,0,0,"B4",Magic+14,0,Green);
      }
//  Si B5 n'existe pas et B1 est en gain de 90 pips, ouvrir B5
      CheckStatus();
      if(BP4!=0 && BP5==0 && Ask>=BP1+90*10*Point) {
      BL5=NL(SL1+SL2+SL3+SL4+SL5+SL6-BL4);
      if(IsTesting()==true) Alert("Hedge B5 = ",BL5,Status);
      BT5=OrderSend(Symbol(),OP_BUY,BL5,Ask,Slip,0,0,"B5",Magic+15,0,Green);
      }
//  Si B6 n'existe pas et B1 est en gain de 110 pips, ouvrir B6
      CheckStatus();
      if(BP5!=0 && BP6==0 && Ask>=BP1+110*10*Point) {
      if(IsTesting()==true) Alert("Hedge B6 = ",BL6,Status);
      BL6=NL(SL1+SL2+SL3+SL4+SL5+SL6-BL5);
      BT6=OrderSend(Symbol(),OP_BUY,BL6,Ask,Slip,0,0,"B6",Magic+16,0,Green);
      }
      }
   }
//+------------------------------------------------------------------+
//| Hedge Sell                                                       |
//+------------------------------------------------------------------+
void HedgeSell()
   {
      if(Bid<=SP1) {
//  Si S2 n'existe pas et S1 est en gain de 30 pips, ouvrir S2
      CheckStatus();
      if(SP1!=0 && SP2==0 && Bid<=SP1-30*10*Point) {
      SL2=NL(SL1);
      if(IsTesting()==true) Alert("Hedge S2 = ",SL2,Status);
      ST2=OrderSend(Symbol(),OP_SELL,SL2,Bid,Slip,0,0,"S2",Magic+22,0,Red);
      }
//  Si S3 n'existe pas et S1 est en gain de 50 pips, ouvrir S3
      CheckStatus();
      if(SP2!=0 && SP3==0 && Bid<=SP1-50*10*Point) {
      SL3=NL(BL1+BL2+BL3+BL4+BL5+BL6-SL2);
      ST3=OrderSend(Symbol(),OP_SELL,SL3,Bid,Slip,0,0,"S3",Magic+23,0,Red);
      if(IsTesting()==true) Alert("Hedge S3 = ",SL3,Status);
      }
//  Si S4 n'existe pas et S1 est en gain de 70 pips, ouvrir S4
      CheckStatus();
      if(SP3!=0 && SP4==0 && Bid<=SP1-70*10*Point) {
      SL4=NL(BL1+BL2+BL3+BL4+BL5+BL6-SL3);
      ST4=OrderSend(Symbol(),OP_SELL,SL4,Bid,Slip,0,0,"S4",Magic+24,0,Red);
      if(IsTesting()==true) Alert("Hedge S4 = ",SL4,Status);
      }
//  Si S5 n'existe pas et S1 est en gain de 90 pips, ouvrir S5
      CheckStatus();
      if(SP4!=0 && SP5==0 && Bid<=SP1-90*10*Point) {
      SL5=NL(BL1+BL2+BL3+BL4+BL5+BL6-SL4);
      ST5=OrderSend(Symbol(),OP_SELL,SL5,Bid,Slip,0,0,"S5",Magic+25,0,Red);
      if(IsTesting()==true) Alert("Hedge S5 = ",SL5,Status);
      }
//  Si S6 n'existe pas et S1 est en gain de 110 pips, ouvrir S6
      CheckStatus();
      if(SP5!=0 && SP6==0 && Bid<=SP1-110*10*Point) {
      SL6=NL(BL1+BL2+BL3+BL4+BL5+BL6-SL5);
      ST6=OrderSend(Symbol(),OP_SELL,SL6,Bid,Slip,0,0,"S6",Magic+26,0,Red);
      if(IsTesting()==true) Alert("Hedge S6 = ",SL6,Status);
      }
      }
   }
//+------------------------------------------------------------------+
//| Set TP Buy                                                       |
//+------------------------------------------------------------------+
void SetTPBuy()
   {
//  Calculer le TP commun à la moyenne des positions plus TP
      CheckStatus();
      if(BP1!=0 && BP2!=0)
      BuyTP=NP(((BP1*BL1+BP2*BL2+BP3*BL3+BP4*BL4+BP5*BL5+BP6*BL6)/(BL1+BL2+BL3+BL4+BL5+BL6))+TP*Pip);
//  Modifier le TP commun des positions ouvertes
      if(OldBuyTP!=BuyTP &&
      MathAbs(BuyTP-Bid)>MinStep && MathAbs(BuyTP-Bid)>Freeze &&
      MathAbs(BP1-Bid)>MinStep && MathAbs(BP1-Bid)>Freeze &&
      MathAbs(BP2-Bid)>MinStep && MathAbs(BP2-Bid)>Freeze &&
      MathAbs(BP3-Bid)>MinStep && MathAbs(BP3-Bid)>Freeze &&
      MathAbs(BP4-Bid)>MinStep && MathAbs(BP4-Bid)>Freeze &&
      MathAbs(BP5-Bid)>MinStep && MathAbs(BP5-Bid)>Freeze &&
      MathAbs(BP6-Bid)>MinStep && MathAbs(BP6-Bid)>Freeze)
      {
      if(BP1!=0 && Bid<BuyTP) Success=OrderModify(BT1,0,0,BuyTP,0,Green);
      if(BP1!=0 && Bid>BuyTP) Success=OrderModify(BT1,0,BuyTP,0,0,Green);
      if(BP2!=0 && Bid<BuyTP) Success=OrderModify(BT2,0,0,BuyTP,0,Green);
      if(BP2!=0 && Bid>BuyTP) Success=OrderModify(BT2,0,BuyTP,0,0,Green);
      if(BP3!=0 && Bid<BuyTP) Success=OrderModify(BT3,0,0,BuyTP,0,Green);
      if(BP3!=0 && Bid>BuyTP) Success=OrderModify(BT3,0,BuyTP,0,0,Green);
      if(BP4!=0 && Bid<BuyTP) Success=OrderModify(BT4,0,0,BuyTP,0,Green);
      if(BP4!=0 && Bid>BuyTP) Success=OrderModify(BT4,0,BuyTP,0,0,Green);
      if(BP5!=0 && Bid<BuyTP) Success=OrderModify(BT5,0,0,BuyTP,0,Green);
      if(BP5!=0 && Bid>BuyTP) Success=OrderModify(BT5,0,BuyTP,0,0,Green);
      if(BP6!=0 && Bid<BuyTP) Success=OrderModify(BT6,0,0,BuyTP,0,Green);
      if(BP6!=0 && Bid>BuyTP) Success=OrderModify(BT6,0,BuyTP,0,0,Green);
//  Enregistrer le TP commun des positions ouvertes
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
//  Calculer le TP commun à la moyenne des positions moins TP
      CheckStatus();
      if(SP1!=0 && SP2!=0)
      SellTP=NP(((SP1*SL1+SP2*SL2+SP3*SL3+SP4*SL4+SP5*SL5+SP6*SL6)/(SL1+SL2+SL3+SL4+SL5+SL6))-TP*Pip);
//  Modifier le TP commun des positions ouvertes
      if(OldSellTP!=SellTP &&
      MathAbs(SellTP-Ask)>MinStep &&
      MathAbs(SellTP-Ask)>Freeze &&
      MathAbs(SP1-Ask)>MinStep &&
      MathAbs(SP2-Ask)>MinStep &&
      MathAbs(SP3-Ask)>MinStep &&
      MathAbs(SP4-Ask)>MinStep &&
      MathAbs(SP5-Ask)>MinStep &&
      MathAbs(SP6-Ask)>MinStep)
      {
      if(SP1!=0 && Ask>SellTP) Success=OrderModify(ST1,0,0,SellTP,0,Red);
      if(SP1!=0 && Ask<SellTP) Success=OrderModify(ST1,0,SellTP,0,0,Red);
      if(SP2!=0 && Ask>SellTP) Success=OrderModify(ST2,0,0,SellTP,0,Red);
      if(SP2!=0 && Ask<SellTP) Success=OrderModify(ST2,0,SellTP,0,0,Red);
      if(SP3!=0 && Ask>SellTP) Success=OrderModify(ST3,0,0,SellTP,0,Red);
      if(SP3!=0 && Ask<SellTP) Success=OrderModify(ST3,0,SellTP,0,0,Red);
      if(SP4!=0 && Ask>SellTP) Success=OrderModify(ST4,0,0,SellTP,0,Red);
      if(SP4!=0 && Ask<SellTP) Success=OrderModify(ST4,0,SellTP,0,0,Red);
      if(SP5!=0 && Ask>SellTP) Success=OrderModify(ST5,0,0,SellTP,0,Red);
      if(SP5!=0 && Ask<SellTP) Success=OrderModify(ST5,0,SellTP,0,0,Red);
      if(SP6!=0 && Ask>SellTP) Success=OrderModify(ST6,0,0,SellTP,0,Red);
      if(SP6!=0 && Ask<SellTP) Success=OrderModify(ST6,0,SellTP,0,0,Red);
//  Enregistrer le TP commun des positions ouvertes
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
//  Si BuyDD < 35% Risk*Cash
      if(BuyDD>-0.35*Risk*Cash) {
//  Si B3 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(BP3!=0 && BG3>=(BL1+BL2+BL3)*10*TP && TimeCurrent()-BO2<=3*24*3600) {
      if(IsTesting()==true) Alert("Turbo B3 = ",BG3);
      Success=OrderClose(BT3,BL3,Bid,Slip,Green);
      }
//  Si B4 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(BP4!=0 && BG4>=(BL1+BL2+BL3+BL4)*10*TP && TimeCurrent()-BO3<=5*24*3600) {
      if(IsTesting()==true) Alert("Turbo B4 = ",BG4);
      Success=OrderClose(BT4,BL4,Bid,Slip,Green);
      }
//  Si B5 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(BP5!=0 && BG5>=(BL1+BL2+BL3+BL4+BL5)*10*TP && TimeCurrent()-BO4<=8*24*3600) {
      if(IsTesting()==true) Alert("Turbo B5 = ",BG5);
      Success=OrderClose(BT5,BL5,Bid,Slip,Green);
      }
      }
   }
//+------------------------------------------------------------------+
//| Turbo Sell                                                       |
//+------------------------------------------------------------------+
void TurboSell()
   {
//  Si SellDD < 35% Risk*Cash
      if(SellDD>-0.35*Risk*Cash) {
//  Si S3 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(SP3!=0 && SG3>=(SL1+SL2+SL3)*10*TP && TimeCurrent()-SO2<=3*24*3600) {
      Success=OrderClose(ST3,SL3,Ask,Slip,Red);
      if(IsTesting()==true) Alert("Turbo S3 = ",SG3);
      }
//  Si S4 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(SP4!=0 && SG4>=(SL1+SL2+SL3+SL4)*10*TP && TimeCurrent()-SO3<=5*24*3600) {
      Success=OrderClose(ST4,SL4,Ask,Slip,Red);
      if(IsTesting()==true) Alert("Turbo S4 = ",SG4);
      }
//  Si S5 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(SP5!=0 && SG5>=(SL1+SL2+SL3+SL4+SL5)*10*TP && TimeCurrent()-SO4<=8*24*3600) {
      Success=OrderClose(ST5,SL5,Ask,Slip,Red);
      if(IsTesting()==true) Alert("Turbo S5 = ",SG5);
      }
      }
   }
//+------------------------------------------------------------------+
//| Close TP Buy                                                     |
//+------------------------------------------------------------------+
void CloseTPBuy()
   {
//  Calculer le TP commun à la moyenne des positions plus TP
      CheckStatus();
      if(BL1!=0)
      BuyTP=NP(((BP1*BL1+BP2*BL2+BP3*BL3+BP4*BL4+BP5*BL5+BP6*BL6)/(BL1+BL2+BL3+BL4+BL5+BL6))+TP*10*Point);
//  Clôturer les positions au TP commun
      if(BL2!=0 && Bid>=BuyTP) {
      if(IsTesting()==true) Alert("Close TP Buy = ",BuyDD,Status);
      CloseBuy();
      }
   }
//+------------------------------------------------------------------+
//| Close TP Sell                                                    |
//+------------------------------------------------------------------+
void CloseTPSell()
   {
//  Calculer le TP commun à la moyenne des positions moins TP
      CheckStatus();
      if(SL1!=0)
      SellTP=NP(((SP1*SL1+SP2*SL2+SP3*SL3+SP4*SL4+SP5*SL5+SP6*SL6)/(SL1+SL2+SL3+SL4+SL5+SL6))-TP*10*Point);
//  Clôturer les positions au TP commun
      if(SL2!=0 && Ask<=SellTP) {
      if(IsTesting()==true) Alert("Close TP Sell = ",SellDD,Status);
      CloseSell();
      }
   }
//+------------------------------------------------------------------+
//| Close Buy                                                        |
//+------------------------------------------------------------------+
void CloseBuy()
   {
//  Clôture des Buy
      CheckStatus();
      if(BP1!=0) Success=OrderClose(BT1,BL1,Bid,Slip,Green);
      if(BP2!=0) Success=OrderClose(BT2,BL2,Bid,Slip,Green);
      if(BP3!=0) Success=OrderClose(BT3,BL3,Bid,Slip,Green);
      if(BP4!=0) Success=OrderClose(BT4,BL4,Bid,Slip,Green);
      if(BP5!=0) Success=OrderClose(BT5,BL5,Bid,Slip,Green);
      if(BP6!=0) Success=OrderClose(BT6,BL6,Bid,Slip,Green);
      if(IsTesting()==true) Alert("CloseBuy = ",BG1+BG2+BG3+BG4+BG5+BG6);
   }
//+------------------------------------------------------------------+
//| Close Sell                                                       |
//+------------------------------------------------------------------+
void CloseSell()
   {
//  Clôture des Sell
      CheckStatus();
      if(SP1!=0) Success=OrderClose(ST1,SL1,Ask,Slip,Red);
      if(SP2!=0) Success=OrderClose(ST2,SL2,Ask,Slip,Red);
      if(SP3!=0) Success=OrderClose(ST3,SL3,Ask,Slip,Red);
      if(SP4!=0) Success=OrderClose(ST4,SL4,Ask,Slip,Red);
      if(SP5!=0) Success=OrderClose(ST5,SL5,Ask,Slip,Red);
      if(SP6!=0) Success=OrderClose(ST6,SL6,Ask,Slip,Red);
      if(IsTesting()==true) Alert("CloseSell = ",SG1+SG2+SG3+SG4+SG5+SG6);
   }
//+------------------------------------------------------------------+
//| Close Buy                                                        |
//+------------------------------------------------------------------+
void CloseBuyOld()
   {
//  Calculer le TP commun à la moyenne des positions ponérées plus TP
      CheckStatus();
      if(BP1!=0)
      BuyTP=NP(((BP1*BL1+BP2*BL2+BP3*BL3+BP4*BL4+BP5*BL5+BP6*BL6)/(BL1+BL2+BL3+BL4+BL5+BL6))+TP*10*Point);
//  Modifier le TP commun des positions ouvertes
      if(Bid>=BuyTP) {
      if(BP1!=0) Success=OrderClose(BT1,BL1,Bid,Slip,Green);
      if(BP2!=0) Success=OrderClose(BT2,BL2,Bid,Slip,Green);
      if(BP3!=0) Success=OrderClose(BT3,BL3,Bid,Slip,Green);
      if(BP4!=0) Success=OrderClose(BT4,BL4,Bid,Slip,Green);
      if(BP5!=0) Success=OrderClose(BT5,BL5,Bid,Slip,Green);
      if(BP6!=0) Success=OrderClose(BT6,BL6,Bid,Slip,Green);
      if(IsTesting()==true) Alert("CloseBuy = ",BG1+BG2+BG3+BG4+BG5+BG6);
      }
   }
//+------------------------------------------------------------------+
//| Close Sell                                                       |
//+------------------------------------------------------------------+
void CloseSellOld()
   {
//  Calculer le TP commun à ma moyenne des positions ponérées moins TP
      CheckStatus();
      if(SP1!=0)
      SellTP=NP(((SP1*SL1+SP2*SL2+SP3*SL3+SP4*SL4+SP5*SL5+SP6*SL6)/(SL1+SL2+SL3+SL4+SL5+SL6))-TP*10*Point);
//  Modifier le TP commun des positions ouvertes
      if(Ask<=SellTP) {
      if(SP1!=0) Success=OrderClose(ST1,SL1,Ask,Slip,Red);
      if(SP2!=0) Success=OrderClose(ST2,SL2,Ask,Slip,Red);
      if(SP3!=0) Success=OrderClose(ST3,SL3,Ask,Slip,Red);
      if(SP4!=0) Success=OrderClose(ST4,SL4,Ask,Slip,Red);
      if(SP5!=0) Success=OrderClose(ST5,SL5,Ask,Slip,Red);
      if(SP6!=0) Success=OrderClose(ST6,SL6,Ask,Slip,Red);
      if(IsTesting()==true) Alert("CloseSell = ",SG1+SG2+SG3+SG4+SG5+SG6);
      }
   }
//+------------------------------------------------------------------+
//| Normalize Price                                                  |
//+------------------------------------------------------------------+
double NP(double price)
   {
//  Normalisation du prix à la cinquième décimale
      price=MathFloor(price*100000)/100000;
      return (price);
   }
//+------------------------------------------------------------------+
//| Normalize Lot                                                    |
//+------------------------------------------------------------------+
double NL(double lot)
   {   
//  Normalisation du lot à la deuxième décimale, lot mini 0.01
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
      
//  Reset des infos sur le graphe
      ObjectsDeleteAll();
//  Position des infos sur le graphe
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

      BarName = "R1";
      ObjectDelete(BarName);
      ObjectCreate(BarName,OBJ_HLINE,0,0,R1);
      ObjectSet(BarName,OBJPROP_COLOR,Yellow);
      ObjectSet(BarName,OBJPROP_WIDTH,1);
      ObjectSet(BarName,OBJPROP_STYLE,3);
      
      BarName = "S1";
      ObjectDelete(BarName);
      ObjectCreate(BarName,OBJ_HLINE,0,0,S1);
      ObjectSet(BarName,OBJPROP_COLOR,Yellow);
      ObjectSet(BarName,OBJPROP_WIDTH,1);
      ObjectSet(BarName,OBJPROP_STYLE,3);
      
      BarName = "R2";
      ObjectDelete(BarName);
      ObjectCreate(BarName,OBJ_HLINE,0,0,R2);
      ObjectSet(BarName,OBJPROP_COLOR,Yellow);
      ObjectSet(BarName,OBJPROP_WIDTH,1);
      ObjectSet(BarName,OBJPROP_STYLE,3);
      
      BarName = "S2";
      ObjectDelete(BarName);
      ObjectCreate(BarName,OBJ_HLINE,0,0,S2);
      ObjectSet(BarName,OBJPROP_COLOR,Yellow);
      ObjectSet(BarName,OBJPROP_WIDTH,1);
      ObjectSet(BarName,OBJPROP_STYLE,3);
      
      BarName = "R3";
      ObjectDelete(BarName);
      ObjectCreate(BarName,OBJ_HLINE,0,0,R3);
      ObjectSet(BarName,OBJPROP_COLOR,Yellow);
      ObjectSet(BarName,OBJPROP_WIDTH,1);
      ObjectSet(BarName,OBJPROP_STYLE,3);
      
      BarName = "S3";
      ObjectDelete(BarName);
      ObjectCreate(BarName,OBJ_HLINE,0,0,S3);
      ObjectSet(BarName,OBJPROP_COLOR,Yellow);
      ObjectSet(BarName,OBJPROP_WIDTH,1);
      ObjectSet(BarName,OBJPROP_STYLE,3);

   }
//+------------------------------------------------------------------+
//| Execute                                                          |
//+------------------------------------------------------------------+
void OnTick()
  {
//  Trader à Minute=0 et Seconde=0 
      if((IsTesting()==true && Seconds()==0) || (IsTesting()==false && IsTradeAllowed()==true)) {
      OpenBuy();
      OpenSell();
//  Trader à Minute=0 et Seconde=0 et si Equity > 50
      if(Trade==true && Minute()==0 && Seconds()==0 && Equity>=50)
//  Si le marché est calme (volatilité M5 < TP pips) depuis 5 minutes
      if(iHigh(Symbol(),PERIOD_M5,1)-iLow(Symbol(),PERIOD_M5,1)<=(TP*10*Point))
//  S'il est moins de 20h
      if(Hour()<20 && Hour()>8) {
      TradeBuy();
      TradeSell();
      }
//  Mettre à jour le TP de sortie des martingales
      SetTPBuy();
      SetTPSell();
      CloseTPBuy();
      CloseTPSell();
//  Couvrir la martingale perdante par une gagnante 
      if(Hedge==true) {
      HedgeBuy();
      HedgeSell();
      }
//  Prendre le gain de la position gagnante
      if(Turbo==true) {
      TurboBuy();
      TurboSell();
      }
//  Clôturer les DD si gain=pertes
      if(Kill==true && Equity>Cash && (BuyDD<-0.45*Risk*Cash || SellDD<-0.45*Risk*Cash)) {
      CloseBuy();
      CloseSell();
      }
//  Afficher les infos des martingales
      Info();
      }
  }
