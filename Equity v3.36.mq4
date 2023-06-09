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
extern double   Risk        =0.3;         // Ratio de risque
extern bool     Start       =true;        // Restart B1 & S1
extern bool     Turbo       =true;        // Profit anticipé
extern bool     Rescue      =true;        // Solder drawdown
extern bool     Unlock      =true;        // Clôture pending
extern bool     Kill        =true;        // Clôture DD

//  Tickets des positions ouvertes
int      BT1, BT2, BT3, BT4, BT5, BT6, BT7;
int      ST1, ST2, ST3, ST4, ST5, ST6, ST7;

//  Prix d'ouverture des positions
double   BP1, BP2, BP3, BP4, BP5, BP6, BP7;
double   SP1, SP2, SP3, SP4, SP5, SP6, SP7;

//  Gain des positions ouvertes
double   BG1, BG2, BG3, BG4, BG5, BG6, BG7;
double   SG1, SG2, SG3, SG4, SG5, SG6, SG7;

//  Volumes des positions ouvertes
double   BL1, BL2, BL3, BL4, BL5, BL6, BL7;
double   SL1, SL2, SL3, SL4, SL5, SL6, SL7;

//  Magic Number des positions ouvertes
int      BM1, BM2, BM3, BM4, BM5, BM6, BM7;
int      SM1, SM2, SM3, SM4, SM5, SM6, SM7;

//  Date d'ouverture des positions ouvertes
datetime BO1, BO2, BO3, BO4, BO5, BO6, BO7;
datetime SO1, SO2, SO3, SO4, SO5, SO6, SO7;

//  TP commun et DD des Buy, Sell
double   BuyTP, SellTP;

//  Divers paramètres
double   Slip=30, Coef=1, Factor=1;
double   Cash=0, Equity=0;
double   BuyDD=0, SellDD=0, BuyDDtime=0, SellDDtime=0;
double   OldBuyTP=0, OldSellTP=0;
double   HW, LW, CW, PP, S1, S2, S3, R1, R2, R3;
bool     Report=true, Success=false;
datetime LastTrade=0;
string   Status=" None";
int      Coin=1;

//+------------------------------------------------------------------+
//| Check Status                                                     |
//+------------------------------------------------------------------+
void CheckStatus()
  {
//  Reset Base de données Buy
      BT1=0; BT2=0; BT3=0; BT4=0; BT5=0; BT6=0; BT7=0;
      BP1=0; BP2=0; BP3=0; BP4=0; BP5=0; BP6=0; BP7=0;
      BG1=0; BG2=0; BG3=0; BG4=0; BG5=0; BG6=0; BG7=0;
      BL1=0; BL2=0; BL3=0; BL4=0; BL5=0; BL6=0; BL7=0;
      BO1=0; BO2=0; BO3=0; BO4=0; BO5=0; BO6=0; BO7=0;

//  Reset Base de données Sell
      ST1=0; ST2=0; ST3=0; ST4=0; ST5=0; ST6=0; ST7=0;
      SP1=0; SP2=0; SP3=0; SP4=0; SP5=0; SP6=0; SP7=0;
      SG1=0; SG2=0; SG3=0; SG4=0; SG5=0; SG6=0; SG7=0;
      SL1=0; SL2=0; SL3=0; SL4=0; SL5=0; SL6=0; SL7=0;
      SO1=0; SO2=0; SO3=0; SO4=0; SO5=0; SO6=0; SO7=0;
      
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
      if(OrderComment()=="B7" || OrderMagicNumber()==Magic+17) {
      BT7=OrderTicket(); BP7=OrderOpenPrice(); BL7=OrderLots(); BG7=MathFloor(OrderProfit()+OrderSwap()); BO7=OrderOpenTime(); }
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
      if(OrderComment()=="S7" || OrderMagicNumber()==Magic+27) {
      ST7=OrderTicket(); SP7=OrderOpenPrice(); SL7=OrderLots(); SG7=MathFloor(OrderProfit()+OrderSwap()); SO7=OrderOpenTime(); }
      }
//  Calcul des DD Buy et Sell
      BuyDD=BG1+BG2+BG3+BG4+BG5+BG6;
      SellDD=SG1+SG2+SG3+SG4+SG5+SG6;
//  Calcul du Cash et de Equity
      Equity=MathRound(AccountEquity()-AccountCredit());
      Cash=MathRound(AccountBalance());
//  Enregistrer le DD des positions ouvertes à chaque heure
      if(Report==true && Minute()==0) {
//  Calcul des points Ouverture, Haut, Bas et Cloture de S et S-1
      HW=MathMax(iHigh(Symbol(),PERIOD_W1,2),iHigh(Symbol(),PERIOD_W1,1));
      LW=MathMin(iLow(Symbol(),PERIOD_W1,2),iLow(Symbol(),PERIOD_W1,1));
      CW=(iClose(Symbol(),PERIOD_W1,2)+iClose(Symbol(),PERIOD_W1,1))/2;
//  Calcul du Point Pivot, Résistances R1, R2, R3 et Supports S1, S2, S3
      PP=(HW+LW+CW)/3;
      R1=2*PP+LW;
      R2=PP+(HW-LW);
      R3=PP+2*(HW-LW);
      S1=2*PP-HW;
      S2=PP-(HW-LW);
      S3=PP-2*(HW-LW);
//  Affichage des infos du compte et des positions
      if(IsTesting()==true)
/*      Alert("Report Price ",Bid," Cash ",Cash," Equity ",Equity,
      " BG1 ",BG1," BG2 ",BG2," BG3 ",BG3," BG4 ",BG4," BG5 ",BG5,
      " BG6 ",BG6," BG7 ",BG7," SG1 ",SG1," SG2 ",SG2," SG3 ",SG3,
      " SG4 ",SG4," SG5 ",SG5," SG6 ",SG6," SG7 ",SG7);
*/
      Report=false;
      }
      if(Minute()!=0) Report=true;
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
//  Si B1 n'existe pas, ouvrir B1
      CheckStatus();
      if(BP1==0) {
      BL1=NL(Factor*0.02);
//  Si BG7 >= 0, BL1 = SLmax/Coef
      if(BG7>=0) {
      if(SL3/Coef>BL1) BL1=NL(SL3/Coef);
      if(SL4/Coef>BL1) BL1=NL(SL4/Coef);
      if(SL5/Coef>BL1) BL1=NL(SL5/Coef);
      }
//  Si BG7 > 10% Cash et 50% BL7 > BL1, BL1 = 50% BL7
      if(BG7>Risk*0.1*Cash && 0.5*BL7>BL1) BL1=NL(0.5*BL7);
      if(IsTesting()==true) Alert("Open B1 = ",BL1,Status);
      BT1=OrderSend(Symbol(),OP_BUY,BL1,Ask,Slip,0,0,"B1",Magic+11,0,Green);
      }
   }
//+------------------------------------------------------------------+
//| Open Sell                                                        |
//+------------------------------------------------------------------+
void OpenSell()
   {
//  Si S1 n'existe pas, ouvrir S1
      CheckStatus();
      if(SP1==0) {
      SL1=NL(Factor*0.02);
//  Si SG7 >= 0, SL1 = BLmax/Coef
      if(SG7>=0) {
      if(BL3/Coef>SL1) SL1=NL(BL3/Coef);
      if(BL4/Coef>SL1) SL1=NL(BL4/Coef);
      if(BL5/Coef>SL1) SL1=NL(BL5/Coef);
      }
//  Si SG7 > 10% Cash et 50% SL7 > SL1, SL1 = 50% SL7
      if(SG7>Risk*0.1*Cash && 0.5*SL7>SL1) SL1=NL(0.5*SL7);
      if(IsTesting()==true) Alert("Open S1 = ",SL1,Status);
      ST1=OrderSend(Symbol(),OP_SELL,SL1,Bid,Slip,0,0,"S1",Magic+21,0,Red);
      }
   }
//+------------------------------------------------------------------+
//| Trade Buy                                                        |
//+------------------------------------------------------------------+
void TradeBuy()
   { 
//  Si BuyDD+SG7 > -35% Cash
      CheckStatus();
      if(BuyDD+SG7>-Risk*0.35*Cash || Unlock==true) {
//  Si B2 n'existe pas et B1 est en perte de 20 pips, ouvrir B2
      CheckStatus();
      if(BP2==0 && BP1!=0 && Ask<=BP1-Coef*20*10*Point) {
      BL2=NL(Factor*0.02);
      if(IsTesting()==true) Alert("Open B2 = ",BL2,Status);
      BT2=OrderSend(Symbol(),OP_BUY,BL2,Ask,Slip,0,0,"B2",Magic+12,0,Green);
      }
//  Si B3 n'existe pas et B2 est en perte de 30 pips, ouvrir B3
      CheckStatus();
      if(BP3==0 && BP2!=0 && Ask<=BP2-Coef*30*10*Point) {
      BL3=NL(Factor*0.03);
      if(IsTesting()==true) Alert("Open B3 = ",BL3,Status);
      BT3=OrderSend(Symbol(),OP_BUY,BL3,Ask,Slip,0,0,"B3",Magic+13,0,Green);
      }
//  Si B4 n'existe pas et B3 est en perte de 50 pips, ouvrir B4
      CheckStatus();
      if(BP4==0 && BP3!=0 && Ask<=BP3-Coef*50*10*Point) {
      BL4=NL(Factor*0.05);
      if(IsTesting()==true) Alert("Open B4 = ",BL4,Status);
      BT4=OrderSend(Symbol(),OP_BUY,BL4,Ask,Slip,0,0,"B4",Magic+14,0,Green);
      }
//  Si B5 n'existe pas et B4 est en perte de 80 pips, ouvrir B5
      CheckStatus();
      if(BP5==0 && BP4!=0 && Ask<=BP4-Coef*80*10*Point) {
      BL5=NL(Factor*0.08);
      if(IsTesting()==true) Alert("Open B5 = ",BL5,Status);
      BT5=OrderSend(Symbol(),OP_BUY,BL5,Ask,Slip,0,0,"B5",Magic+15,0,Green);
      }
//  Si B6 n'existe pas et B5 est en perte de 120 pips, ouvrir B6
      CheckStatus();
      if(BP6==0 && BP5!=0 && Ask<=BP5-Coef*130*10*Point) {
      if(IsTesting()==true) Alert("Open B6 = ",BL6,Status);
      BL6=NL(Factor*0.13);
      BT6=OrderSend(Symbol(),OP_BUY,BL6,Ask,Slip,0,0,"B6",Magic+16,0,Green);
      }
      }
   }
//+------------------------------------------------------------------+
//| Trade Sell                                                       |
//+------------------------------------------------------------------+
void TradeSell()
   {
//  Si SellDD+BG7 > -35% Cash
      CheckStatus();
      if(SellDD+BG7>-Risk*0.35*Cash || Unlock==true) {
//  Si S2 n'existe pas et S1 est en perte de 20 pips, ouvrir S2
      CheckStatus();
      if(SP2==0 && SP1!=0 && Bid>=SP1+Coef*20*10*Point) {
      SL2=NL(Factor*0.02);
      if(IsTesting()==true) Alert("Open S2 = ",SL2,Status);
      ST2=OrderSend(Symbol(),OP_SELL,SL2,Bid,Slip,0,0,"S2",Magic+22,0,Red);
      }
//  Si S3 n'existe pas et S2 est en perte de 30 pips, ouvrir S3
      CheckStatus();
      if(SP3==0 && SP2!=0 && Bid>=SP2+Coef*30*10*Point) {
      SL3=NL(Factor*0.03);
      if(IsTesting()==true) Alert("Open S3 = ",SL3,Status);
      ST3=OrderSend(Symbol(),OP_SELL,SL3,Bid,Slip,0,0,"S3",Magic+23,0,Red);
      }
//  Si S4 n'existe pas et S3 est en perte de 50 pips, ouvrir S4
      CheckStatus();
      if(SP4==0 && SP3!=0 && Bid>=SP3+Coef*50*10*Point) {
      SL4=NL(Factor*0.05);
      if(IsTesting()==true) Alert("Open S4 = ",SL4,Status);
      ST4=OrderSend(Symbol(),OP_SELL,SL4,Bid,Slip,0,0,"S4",Magic+24,0,Red);
      }
//  Si S5 n'existe pas et S4 est en perte de 80 pips, ouvrir S5
      CheckStatus();
      if(SP5==0 && SP4!=0 && Bid>=SP4+Coef*80*10*Point) {
      SL5=NL(Factor*0.08);
      if(IsTesting()==true) Alert("Open S5 = ",SL5,Status);
      ST5=OrderSend(Symbol(),OP_SELL,SL5,Bid,Slip,0,0,"S5",Magic+25,0,Red);
      }
//  Si S6 n'existe pas et S5 est en perte de 120 pips, ouvrir S6
      CheckStatus();
      if(SP6==0 && SP5!=0 && Bid>=SP4+Coef*130*10*Point) {
      SL6=NL(Factor*0.13);
      if(IsTesting()==true) Alert("Open S6 = ",SL6,Status);
      ST6=OrderSend(Symbol(),OP_SELL,SL6,Bid,Slip,0,0,"S6",Magic+26,0,Red);
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
      if(BL1!=0) {
      BuyTP=NP(((BP1*BL1+BP2*BL2+BP3*BL3+BP4*BL4+BP5*BL5+BP6*BL6+BP7*BL7)/(BL1+BL2+BL3+BL4+BL5+BL6+BL7))+TP*10*Point);
//  Si B7 existe et BP7 < BuyTP, recalculer le TP commun sans B7
      if(BP7!=0 && BP7<BuyTP)
      BuyTP=NP(((BP1*BL1+BP2*BL2+BP3*BL3+BP4*BL4+BP5*BL5+BP6*BL6)/(BL1+BL2+BL3+BL4+BL5+BL6))+TP*10*Point);
      }
//  Modifier le TP commun des positions ouvertes
      if(OldBuyTP!=BuyTP) { 
      if(BP1!=0) Success=OrderModify(BT1,0,0,BuyTP,0,Green);
      if(BP2!=0) Success=OrderModify(BT2,0,0,BuyTP,0,Green);
      if(BP3!=0) Success=OrderModify(BT3,0,0,BuyTP,0,Green);
      if(BP4!=0) Success=OrderModify(BT4,0,0,BuyTP,0,Green);
      if(BP5!=0) Success=OrderModify(BT5,0,0,BuyTP,0,Green);
      if(BP6!=0) Success=OrderModify(BT6,0,0,BuyTP,0,Green);
//  Si B7 existe et BuyTP < BP7, modifier aussi le TP B7
      if(BG7!=0 && BP7>BuyTP) Success=OrderModify(BT7,0,0,BuyTP,0,Green);
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
      if(SL1!=0) {
      SellTP=NP(((SP1*SL1+SP2*SL2+SP3*SL3+SP4*SL4+SP5*SL5+SP6*SL6+SP7*SL7)/(SL1+SL2+SL3+SL4+SL5+SL6+SL7))-TP*10*Point);
//  Si S7 existe et SP7 > SellTP, recalculer le TP commun sans S7
      if(SP7!=0 && SP7>SellTP)
      SellTP=NP(((SP1*SL1+SP2*SL2+SP3*SL3+SP4*SL4+SP5*SL5+SP6*SL6)/(SL1+SL2+SL3+SL4+SL5+SL6))-TP*10*Point);
      }
//  Modifier le TP commun des positions ouvertes
      if(OldSellTP!=SellTP) {
      if(SP1!=0) Success=OrderModify(ST1,0,0,SellTP,0,Red);
      if(SP2!=0) Success=OrderModify(ST2,0,0,SellTP,0,Red);
      if(SP3!=0) Success=OrderModify(ST3,0,0,SellTP,0,Red);
      if(SP4!=0) Success=OrderModify(ST4,0,0,SellTP,0,Red);
      if(SP5!=0) Success=OrderModify(ST5,0,0,SellTP,0,Red);
      if(SP6!=0) Success=OrderModify(ST6,0,0,SellTP,0,Red);
//  Si S7 existe et SellTP > SP7, modifier aussi le TP S7
      if(SP7!=0 && SP7<SellTP) Success=OrderModify(ST7,0,0,SellTP,0,Red);
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
//  Si SG7 >= 0
      CheckStatus();
      if(SG7>=0) {
//  Si B2 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(BP2!=0 && BG2>=(BL1+BL2)*10*TP && TimeCurrent()-BO1<=2*24*3600) {
      if(IsTesting()==true) Alert("Turbo B2 = ",BG2);
      Success=OrderClose(BT2,BL2,Bid,Slip,Green);
      }
//  Si B3 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(BP3!=0 && BG3>=(BL1+BL2+BL3)*10*TP && TimeCurrent()-BO2<=2*24*3600) {
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
//  Si BG7 >= 0
      CheckStatus();
      if(BG7>=0) {
//  Si S2 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(SP2!=0 && SG2>=(SL1+SL2)*10*TP && TimeCurrent()-SO1<=2*24*3600) {
      if(IsTesting()==true) Alert("Turbo S2 = ",SG2);
      Success=OrderClose(ST2,SL2,Ask,Slip,Red);
      }
//  Si S3 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(SP3!=0 && SG3>=(SL1+SL2+SL3)*10*TP && TimeCurrent()-SO2<=2*24*3600) {
      if(IsTesting()==true) Alert("Turbo S3 = ",SG3);
      Success=OrderClose(ST3,SL3,Ask,Slip,Red);
      }
//  Si S4 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(SP4!=0 && SG4>=(SL1+SL2+SL3+SL4)*10*TP && TimeCurrent()-SO3<=5*24*3600) {
      if(IsTesting()==true) Alert("Turbo S4 = ",SG4);
      Success=OrderClose(ST4,SL4,Ask,Slip,Red);
      }
//  Si S5 fait le gain du TP commun, la cloturer
      CheckStatus();
      if(SP5!=0 && SG5>=(SL1+SL2+SL3+SL4+SL5)*10*TP && TimeCurrent()-SO4<=8*24*3600) {
      if(IsTesting()==true) Alert("Turbo S5 = ",SG5);
      Success=OrderClose(ST5,SL5,Ask,Slip,Red);
      }
      }
   }
//+------------------------------------------------------------------+
//| Hedge Buy                                                        |
//+------------------------------------------------------------------+
void HedgeBuy()
   {
//  Si S7 n'existe pas et BuyDD < -35% Cash, ouvrir S7
      CheckStatus();
      if(SP7==0 && BuyDD<=-Risk*0.35*Cash) {
//  Calculer SL7 = 80% (BL1+BL2+BL3+BL4+BL5+BL6+BL7)
      SL7=NL(0.8*(BL1+BL2+BL3+BL4+BL5+BL6+BL7));
      if(IsTesting()==true) Alert("Open S7 = ",SL7,Status);
      ST7=OrderSend(Symbol(),OP_SELL,SL7,Bid,Slip,0,0,"S7",Magic+27,0,Yellow);
      BuyDDtime=TimeCurrent();
      }
   }
//+------------------------------------------------------------------+
//| Hedge Sell                                                       |
//+------------------------------------------------------------------+
void HedgeSell()
   {
//  Si B7 n'existe pas et SellDD < -35% Cash, ouvrir B7
      CheckStatus();
      if(BP7==0 && SellDD<=-Risk*0.35*Cash) {
//  Calculer BL7 = 80% (SL1+SL2+SL3+SL4+SL5+SL6+SL7)
      BL7=NL(0.8*(SL1+SL2+SL3+SL4+SL5+SL6+SL7));
      if(IsTesting()==true) Alert("Open B7 = ",BL7,Status);
      BT7=OrderSend(Symbol(),OP_BUY,BL7,Ask,Slip,0,0,"B7",Magic+17,0,Yellow);
      SellDDtime=TimeCurrent();
      }
   }
//+------------------------------------------------------------------+
//| Rescue Buy                                                       |
//+------------------------------------------------------------------+
void RescueBuy()
   {
//  Si B6 et B7 existent et SG1+BG6+BG7+SG7 > (SL1+BL6+BL7+SL7)*10*TP, les cloturer
      CheckStatus();
      if(BP6!=0 && BP7!=0 && SP7!=0 && SG1+BG6+BG7+SG7>(SL1+BL6+BL7+SL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue S1+B6+B7+S7 = ",SG1+BG6+BG7+SG7,Status);
      Success=OrderClose(ST1,SL1,Ask,Slip,Red);
      Success=OrderClose(BT6,BL6,Bid,Slip,Green);
      Success=OrderClose(BT7,BL7,Bid,Slip,Green);
      Success=OrderClose(ST7,SL7,Ask,Slip,Red);
      }
//  Si S7 existe mais pas B7 et
      CheckStatus();
      if(SP7!=0 && BP7==0) {
//  Si B5 et B6 existent mais pas B7 et SG1+BG5+BG6+SG7 > (SL1+BL5+BL6+SL7)*10*TP, les cloturer
      CheckStatus();
      if(BP5!=0 && BP6!=0 && SG1+BG5+BG6+SG7>(SL1+BL5+BL6+SL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue S1+B5+B6+S7 = ",SG1+BG5+BG6+SG7,Status);
      Success=OrderClose(ST1,SL1,Ask,Slip,Red);
      Success=OrderClose(BT5,BL5,Bid,Slip,Green);
      Success=OrderClose(BT6,BL6,Bid,Slip,Green);
      Success=OrderClose(ST7,SL7,Ask,Slip,Red);
      }
//  Si B4 et B6 existent mais pas B5 et SG1+BG4+BG6+SG7 > (SL1+BL4+BL6+SL7)*10*TP, les cloturer
      CheckStatus();
      if(BP4!=0 && BP5==0 && BP6!=0 && SG1+BG4+BG6+SG7>(SL1+BL4+BL6+SL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue S1+B4+B6+S7 = ",SG1+BG4+BG6+SG7,Status);
      Success=OrderClose(ST1,SL1,Ask,Slip,Red);
      Success=OrderClose(BT4,BL4,Bid,Slip,Green);
      Success=OrderClose(BT6,BL6,Bid,Slip,Green);
      Success=OrderClose(ST7,SL7,Ask,Slip,Red);
      }
//  Si B3 et B6 existent mais pas B4 et SG1+BG3+BG6+SG7 > (SL1+BL3+BL6+SL7)*10*TP, les cloturer
      CheckStatus();
      if(BP3!=0 && BP4==0 && BP6!=0 && SG1+BG3+BG6+SG7>(SL1+BL3+BL6+SL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue S1+B3+B6+S7 = ",SG1+BG3+BG6+SG7,Status);
      Success=OrderClose(ST1,SL1,Ask,Slip,Red);
      Success=OrderClose(BT3,BL3,Bid,Slip,Green);
      Success=OrderClose(BT6,BL6,Bid,Slip,Green);
      Success=OrderClose(ST7,SL7,Ask,Slip,Red);
      }
//  Si B6 existe mais pas B3 et SG1+BG3+BG6+SG7 > (SL1+BL6+SL7)*10*TP, les cloturer
      CheckStatus();
      if(BP3==0 && BP6!=0 && SG1+BG6+SG7>(SL1+BL6+SL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue S1+B6+S7 = ",SG1+BG6+SG7,Status);
      Success=OrderClose(ST1,SL1,Ask,Slip,Red);
      Success=OrderClose(BT6,BL6,Bid,Slip,Green);
      Success=OrderClose(ST7,SL7,Ask,Slip,Red);
      }
      }
   }
//+------------------------------------------------------------------+
//| Rescue Sell                                                      |
//+------------------------------------------------------------------+
void RescueSell()
   {
//  Si S6, S7 et B7 existent et BG1+SG6+SG7+BG7 > (BL1+SL6+SL7+BL7)*10*TP, les clôturer
      CheckStatus();
      if(SP6!=0 && SP7!=0 && BP7!=0 && BG1+SG6+SG7+BG7>(BL1+SL6+SL7+BL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue B1+S6+S7+B7 = ",BG1+SG6+SG7+BG7,Status);
      Success=OrderClose(BT1,BL1,Bid,Slip,Green);
      Success=OrderClose(ST6,SL6,Ask,Slip,Red);
      Success=OrderClose(ST7,SL7,Ask,Slip,Red);
      Success=OrderClose(BT7,BL7,Bid,Slip,Green);
      }
//  Si B7 existe mais pas S7 et
      CheckStatus();
      if(BP7!=0 && SP7==0) {
//  Si S5 et S6 existent et BG1+SG5+SG6+BG7 > (BL1+SL5+SL6+BL7)*10*TP, les clôturer
      CheckStatus();
      if(SP5!=0 && SP6!=0 && BG1+SG5+SG6+BG7>(BL1+SL5+SL6+BL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue B1+S5+S6+B7 = ",BG1+SG5+SG6+BG7,Status);
      Success=OrderClose(BT1,BL1,Bid,Slip,Green);
      Success=OrderClose(ST5,SL5,Ask,Slip,Red);
      Success=OrderClose(ST6,SL6,Ask,Slip,Red);
      Success=OrderClose(BT7,BL7,Bid,Slip,Green);
      }
//  Si S4 et S6 existent mais pas S5 et BG1+SG4+SG6+BG7 > (BL1+SL4+SL6+BL7)*10*TP, les clôturer
      CheckStatus();
      if(SP4!=0 && SP6!=0 && SP5==0 && BG1+SG4+SG6+BG7>(BL1+SL4+SL6+BL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue B1+S4+S6+B7 = ",BG1+SG4+SG6+BG7,Status);
      Success=OrderClose(BT1,BL1,Bid,Slip,Green);
      Success=OrderClose(ST4,SL4,Ask,Slip,Red);
      Success=OrderClose(ST6,SL6,Ask,Slip,Red);
      Success=OrderClose(BT7,BL7,Bid,Slip,Green);
      }
//  Si S3 et S6 existent mais pas S4 et BG1+SG3+SG6+BG7 > (BL1+SL3+SL6+BL7)*10*TP, les clôturer
      CheckStatus();
      if(SP3!=0 && SP6!=0 && SP4==0 && BG1+SG3+SG6+BG7>(BL1+SL3+SL6+BL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue B1+S3+S6+B7 = ",BG1+SG3+SG6+BG7,Status);
      Success=OrderClose(BT1,BL1,Bid,Slip,Green);
      Success=OrderClose(ST3,SL3,Ask,Slip,Red);
      Success=OrderClose(ST6,SL6,Ask,Slip,Red);
      Success=OrderClose(BT7,BL7,Bid,Slip,Green);
      }
//  Si S6 existe mais pas S3 et BG1+SG6+BG7 > (BL1+SL6+BL7)*10*TP, les clôturer
      CheckStatus();
      if(SP6!=0 && SP3==0 && BG1+SG3+SG6+BG7>(BL1+SL6+BL7)*10*TP) {
      if(IsTesting()==true) Alert("Rescue B1+S6+B7 = ",BG1+SG6+BG7,Status);
      Success=OrderClose(BT1,BL1,Bid,Slip,Green);
      Success=OrderClose(ST6,SL6,Ask,Slip,Red);
      Success=OrderClose(BT7,BL7,Bid,Slip,Green);
      }
      }
   }
//+------------------------------------------------------------------+
//| Unlock Buy                                                       |
//+------------------------------------------------------------------+
void UnlockBuy()
   {
//  Si B2 et B3 existent et BG2+BG3 > (BL2+BL3)*10*TP et B2 > 1j, les cloturer
      CheckStatus();
      if(BP2!=0 && BP3!=0 && BG2+BG3>(BL2+BL3)*10*TP && TimeCurrent()-BO2>=1*24*3600) {
      if(IsTesting()==true) Alert("Rescue B2+B3 = ",BG2+BG3,Status);
      Success=OrderClose(BT2,BL2,Bid,Slip,Green);
      Success=OrderClose(BT3,BL3,Bid,Slip,Green);
      }
//  Si B3 et B4 existent et BG3+BG4 > (BL3+BL4)*10*TP et B3 > 2j, les cloturer
      CheckStatus();
      if(BP3!=0 && BP4!=0 && BG3+BG4>(BL3+BL4)*10*TP && TimeCurrent()-BO3>=2*24*3600) {
      if(IsTesting()==true) Alert("Rescue B3+B4 = ",BG3+BG4,Status);
      Success=OrderClose(BT3,BL3,Bid,Slip,Green);
      Success=OrderClose(BT4,BL4,Bid,Slip,Green);
      }
//  Si B4 et B5 existent et BG4+BG5 > (BL4+BL5)*10*TP et B42 > 5j, les cloturer
      CheckStatus();
      if(BP4!=0 && BP5!=0 && BG4+BG5>(BL4+BL5)*10*TP && TimeCurrent()-BO4>=5*24*3600) {
      if(IsTesting()==true) Alert("Rescue B4+B5 = ",BG4+BG5,Status);
      Success=OrderClose(BT4,BL4,Bid,Slip,Green);
      Success=OrderClose(BT5,BL5,Bid,Slip,Green);
      }
//  Si B5 et B6 existent et BG5+BG6 > (BL5+BL6)*10*TP et B52 > 8j, les cloturer
      CheckStatus();
      if(BP5!=0 && BP6!=0 && BG5+BG6>(BL5+BL6)*10*TP && TimeCurrent()-BO5>=8*24*3600) {
      if(IsTesting()==true) Alert("Rescue B5+B6 = ",BG5+BG6,Status);
      Success=OrderClose(BT5,BL5,Bid,Slip,Green);
      Success=OrderClose(BT6,BL6,Bid,Slip,Green);
      }
   }
//+------------------------------------------------------------------+
//| Unlock Sell                                                      |
//+------------------------------------------------------------------+
void UnlockSell()
   {
//  Si S2 et S3 existent et SG2+SG3 > (SL2+SL3)*10*TP et SB2 > 1j, les cloturer
      CheckStatus();
      if(SP2!=0 && SP3!=0 && SG2+SG3>(SL2+SL3)*10*TP && TimeCurrent()-SO2>=1*24*3600) {
      if(IsTesting()==true) Alert("Rescue S2+S3 = ",SG2+SG3,Status);
      Success=OrderClose(ST2,SL2,Ask,Slip,Red);
      Success=OrderClose(ST3,SL3,Ask,Slip,Red);
      }
//  Si S3 et S4 existent et SG3+SG4 > (SL3+SL4)*10*TP et SB2 > 2j, les cloturer
      CheckStatus();
      if(SP3!=0 && SP4!=0 && SG3+SG4>(SL3+SL4)*10*TP && TimeCurrent()-SO3>=2*24*3600) {
      if(IsTesting()==true) Alert("Rescue S3+S4 = ",SG3+SG4,Status);
      Success=OrderClose(ST3,SL3,Ask,Slip,Red);
      Success=OrderClose(ST4,SL4,Ask,Slip,Red);
      }
//  Si S4 et S5 existent et SG4+SG5 > (SL4+SL5)*10*TP et SB2 > 5j, les cloturer
      CheckStatus();
      if(SP4!=0 && SP5!=0 && SG4+SG5>(SL4+BL5)*10*TP && TimeCurrent()-SO4>=5*24*3600) {
      if(IsTesting()==true) Alert("Rescue S4+S5 = ",SG4+SG5,Status);
      Success=OrderClose(ST4,SL4,Ask,Slip,Red);
      Success=OrderClose(ST5,SL5,Ask,Slip,Red);
      }
//  Si S5 et S6 existent et SG5+SG6 > (SL5+SL6)*10*TP et SB2 > 8j, les cloturer
      CheckStatus();
      if(SP5!=0 && SP6!=0 && SG5+SG6>(SL5+SL6)*10*TP && TimeCurrent()-SO5>=8*24*3600) {
      if(IsTesting()==true) Alert("Rescue S5+S6 = ",SG5+SG6,Status);
      Success=OrderClose(ST5,SL5,Ask,Slip,Red);
      Success=OrderClose(ST6,SL6,Ask,Slip,Red);
      }
   }
//+------------------------------------------------------------------+
//| Kill BuyDD                                                       |
//+------------------------------------------------------------------+
void KillBuyDD()
   {
//  Si SellDD < BuyDD et SellDDtime > 30 jours, cloturer SellDD
      CheckStatus();
      if(BuyDD < SellDD && TimeCurrent()-BuyDDtime>=30*24*3600) {
      if(SP7!=0) Success=OrderClose(ST7,SL7,Ask,Slip,Red);
      CloseBuy();
      }
   }
//+------------------------------------------------------------------+
//| Kill SellDD                                                      |
//+------------------------------------------------------------------+
void KillSellDD()
   {
//  Si BuyDD < SellDD et BuyDDtime > 30 jours, cloturer BuyDD
      CheckStatus();
      if(SellDD < BuyDD && TimeCurrent()-SellDDtime>=30*24*3600) {
      if(BP7!=0) Success=OrderClose(BT7,BL7,Bid,Slip,Green);
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
      ObjectSetText(BarName,"Risk="+NL(Risk)+" Coef="+NL(Coef)+" Factor="+NL(Factor)+" TP="+TP,12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,20);
      ObjectSet(BarName,OBJPROP_YDISTANCE,20);

      CheckStatus();
      BarName = "Parameters2";
      ObjectCreate(BarName,OBJ_LABEL,0,0,0);
      ObjectSetText(BarName,"Turbo="+Turbo+" Rescue="+Rescue+" Unlock="+Unlock+" Kill="+Kill,12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,20);
      ObjectSet(BarName,OBJPROP_YDISTANCE,40);

      BarName = "DD";
      ObjectCreate(BarName,OBJ_LABEL,0,0,0);
      ObjectSetText(BarName,"BDD= "+MathRound(100*(BuyDD+BG7)/Cash)+" %"+" SDD = "+MathRound(100*(SellDD+SG7)/Cash)+" %",12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,20);
      ObjectSet(BarName,OBJPROP_YDISTANCE,60);

      CheckStatus();
      BarName = "Money";
      ObjectCreate(BarName,OBJ_LABEL,0,0,0);
      ObjectSetText(BarName,"Cash= "+MathRound(Cash)+"€"+" Equity= "+MathRound(AccountEquity())+"€",12,"Corbel",YellowGreen);
      ObjectSet(BarName,OBJPROP_CORNER,Coin);
      ObjectSet(BarName,OBJPROP_XDISTANCE,20);
      ObjectSet(BarName,OBJPROP_YDISTANCE,80);

      BarName = "PP";
      ObjectDelete(BarName);
      ObjectCreate(BarName,OBJ_HLINE,0,0,PP);
      ObjectSet(BarName,OBJPROP_COLOR,Blue);
      ObjectSet(BarName,OBJPROP_WIDTH,2);
      ObjectSet(BarName,OBJPROP_STYLE,0);

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
//  Ouvrir les positions initiales
      if(Start==true) {
      OpenBuy();
      OpenSell();
      }
//  Mettre à jour le TP de sortie des martingales
      SetTPBuy();
      SetTPSell();
//  Couvrir les martingales perdantes
      HedgeBuy();
      HedgeSell();
//  Récupérer les positions perdantes
      if(Rescue==true) {
      RescueBuy();
      RescueSell();
      }
//  Récupérer les positions perdantes
      if(Unlock==true) {
      UnlockBuy();
      UnlockSell();
      }

//  Trader à Minute=0 et Seconde=0 et si Equity > 250
      if(Minute()==0 && Seconds()==0 && Equity>=35)
//  Si le marché est calme (volatilité M5 < 10 TP) depuis 5 minutes
      if(iHigh(Symbol(),PERIOD_M5,1)-iLow(Symbol(),PERIOD_M5,1)<=(TP*10*Point))
//  S'il est moins de 20h
      if(Hour()<20 && Hour()>8) {
      TradeBuy();
      TradeSell();
      }
//  Clôturer les positions gagnantes
      if(Turbo==true) {
      TurboBuy();
      TurboSell();
      }
//  Afficher les infos des martingales
      Info();
      
      }
  }
