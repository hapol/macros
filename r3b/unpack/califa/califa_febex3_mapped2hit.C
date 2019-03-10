/* Additional info:
 * To generate the header file for a different number of febex modules,
 * add then in $UCESB_DIR/upexps/califalisbon16/califa.spec file and use:
 *
 * $UCESB_DIR/upexps/califalisbon16/califa --ntuple=RAW,STRUCT_HH,ext_h101_raw_califa_febex.h
 *
 * Put this header file into the 'r3bsource' directory and recompile
 * (present file is valid for 9 modules and 16 channels for module)
 * */

 typedef struct EXT_STR_h101_t {
   EXT_STR_h101_unpack_t unpack;
   EXT_STR_h101_CALIFA_t califa;
 } EXT_STR_h101;

void califa_febex3_mapped2cal() {
  TStopwatch timer;
  timer.Start();

  //const Int_t nev = -1; /* number of events to read, -1 - until CTRL+C */
  const Int_t nev = -1; /* number of events to read, -1 - until CTRL+C */

  /* Create source using ucesb for input ------------------ */
  TString filename = "/Users/hapol/data_0002.lmd";
  TString outputFileName = "/Users/hapol/hits_0002.root";
  TString ntuple_options = "UNPACK:EVENTNO,UNPACK:TRIGGER,RAW";
  TString ucesb_dir = getenv("UCESB_DIR");
  TString ucesb_path = ucesb_dir + "/../upexps/201902_s444";
  EXT_STR_h101 ucesb_struct;

  R3BUcesbSource* source = new R3BUcesbSource(filename, ntuple_options,
					      ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
  source->SetMaxEvents(nev);

  source->AddReader(new R3BUnpackReader((EXT_STR_h101_unpack*)&ucesb_struct,
					offsetof(EXT_STR_h101, unpack)));
  source->AddReader(new R3BCalifaFebexReader((EXT_STR_h101_CALIFA*)&ucesb_struct.califa,
					     offsetof(EXT_STR_h101, califa)));
  /* ------------------------------------------------------ */

  /* Create online run ------------------------------------ */
  FairRunOnline* run = new FairRunOnline(source);
  run->SetRunId(1513078509);
  run->SetOutputFile(outputFileName);

  /* Add analysis task ------------------------------------ */

  //R3BCalifaMapped2CrystalCal ---
  R3BCalifaMapped2CrystalCal* Map2Cal = new R3BCalifaMapped2CrystalCal();
  run->AddTask(Map2Cal);

  //R3BCalifaCrystalCal2Hit ---
  R3BCalifaCrystalCal2Hit* califaHF = new R3BCalifaCrystalCal2Hit();
  TMap fDetList = new TMap();
  califaHF->SetGeometryFileName("/Users/hapol/CODES/R3BROOT/R3BRoot/geometry/califa_s444.geo.root");
  califaHF->SelectGeometryVersion(444);
  califaHF->SetClusteringAlgorithm(1,0);
  califaHF->SetDetectionThreshold(0.000050);  //50 KeV
  califaHF->SetExperimentalResolution(0.);    //6% at 1 MeV
  //califaHF->SetComponentResolution(0.);      //sigma = 0.25 MeV
  //califaHF->SetPhoswichResolution(3.,5.);     //percent @ 1 MeV for LaBr and LaCl
  califaHF->SetAngularWindow(0.2,0.2);        //[0.25 around 14.3 degrees, 3.2 for the complete calorimeter]
  run->AddTask(califaHF);


  /* Initialize ------------------------------------------- */
  //run->Init();
  FairLogger::GetLogger()->SetLogScreenLevel("INFO");
  /* ------------------------------------------------------ */

  /* Runtime data base ------------------------------------ */
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  //FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();//Ascii
  //parIo1->open("Params.par","in");
  FairParRootFileIo* parIo1 = new FairParRootFileIo();//Root
  parIo1->open("Califa_CalibParam.root","in");
  rtdb->setFirstInput(parIo1);
  rtdb->print();
  /* ------------------------------------------------------ */

  /* Run -------------------------------------------------- */
  run->Init();
  run->Run(nev,0);
  //rtdb->saveOutput();
  delete run;
  /* ------------------------------------------------------ */

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Output file is " << outputFileName << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s"
       << endl << endl;
}
