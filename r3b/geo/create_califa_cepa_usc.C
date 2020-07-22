// CALIFA CEPA, USC VERSION (JUL 2020)
// Requires files CLF-simu-CEPA_USC.txt
#include "TGeoManager.h"
#include "TMath.h"
#include <iomanip>
#include <iostream>

// Create Matrix Unity
TGeoRotation* fGlobalRot = new TGeoRotation();

// Create a null translation
TGeoTranslation* fGlobalTrans = new TGeoTranslation();
TGeoRotation* fRefRot = NULL;

TGeoManager* gGeoMan = NULL;

Double_t fThetaX = 0.;
Double_t fThetaY = 0.;
Double_t fThetaZ = 0.;
Double_t fPhi = 0.;
Double_t fTheta = 0.;
Double_t fPsi = 0.;
Double_t fX = 0.;
Double_t fY = 0.;
Double_t fZ = 0.;
Bool_t fLocalTrans = kFALSE;
Bool_t fLabTrans = kFALSE;

TGeoCombiTrans* GetGlobalPosition(TGeoCombiTrans* fRef);

Bool_t isCrystalInstalled(Int_t alvType, Int_t alveolusCopy);

void create_califa_geo(const char* geoTag = "")
{
    fGlobalTrans->SetTranslation(0.0, 0.0, 0.0);

    // -------   Load media from media file   -------------------------
    FairGeoLoader* geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
    FairGeoInterface* geoFace = geoLoad->getGeoInterface();
    TString geoPath = gSystem->Getenv("VMCWORKDIR");
    TString medFile = geoPath + "/geometry/media_r3b.geo";
    geoFace->setMediaFile(medFile);
    geoFace->readMedia();
    gGeoMan = gGeoManager;

    // -------   Geometry file name (output)   ------------------------
    TString geoFileName = geoPath + "/geometry/califa_";
    geoFileName = geoFileName + geoTag + ".geo.root";

    // -----------------   Get and create the required media    -------
    FairGeoMedia* geoMedia = geoFace->getMedia();
    FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

    FairGeoMedium* mAir = geoMedia->getMedium("Air");
    if (!mAir)
        Fatal("Main", "FairMedium Air not found");
    geoBuild->createMedium(mAir);
    TGeoMedium* pAirMedium = gGeoMan->GetMedium("Air");
    if (!pAirMedium)
        Fatal("Main", "Medium Air not found");

    FairGeoMedium* mCsI = geoMedia->getMedium("CsI");
    if (!mCsI)
        Fatal("Main", "FairMedium CsI not found");
    geoBuild->createMedium(mCsI);
    TGeoMedium* pCsIMedium = gGeoMan->GetMedium("CsI");
    if (!pCsIMedium)
        Fatal("Main", "Medium CsI not found");

    FairGeoMedium* mCar = geoMedia->getMedium("CarbonFibre");
    if (!mCar)
        Fatal("Main", "FairMedium CarbonFibre not found");
    geoBuild->createMedium(mCar);
    TGeoMedium* pCarbonFibreMedium = gGeoMan->GetMedium("CarbonFibre");
    if (!pCarbonFibreMedium)
        Fatal("Main", "Medium CarbonFibre not found");

    FairGeoMedium* mTfl = geoMedia->getMedium("Tefflon");
    if (!mTfl)
        Fatal("Main", "FairMedium Tefflon not found");
    geoBuild->createMedium(mTfl);
    TGeoMedium* pWrappingMedium = gGeoMan->GetMedium("Tefflon");
    if (!pWrappingMedium)
        Fatal("Main", "Medium Tefflon not found");
    // ----------------------------------------------------------------

    // --------------   Create geometry and top volume  ---------------
    gGeoMan = (TGeoManager*)gROOT->FindObject("FAIRGeom");
    gGeoMan->SetName("CALIFAgeom");
    TGeoVolume* top = new TGeoVolumeAssembly("TOP");
    gGeoMan->SetTopVolume(top);
    // ----------------------------------------------------------------

    // Defintion of the Mother Volume
    Double_t length = 300.;
    TGeoShape* pCBWorldOut = new TGeoTube("Califa_boxOut",
                                          0.,    // Rmin
                                          100.,  // Rmax
                                          100.); // half length

    TGeoShape* pCBWorldIn1 = new TGeoTube("Califa_Centerpart1", // hole to accommodate the tracker
                                          0.,                   // Rmin
                                          26.4,                 // Rmax
                                          65.);                 // half length
    TGeoCombiTrans* t_part1 = new TGeoCombiTrans("t_part1", 0., 0., -35, fRefRot);
    t_part1->RegisterYourself();

    TGeoShape* pCBWorldIn2 = new TGeoTube("Califa_Centerpart2", // hole to accommodate the pipe through the end-cap
                                          0.,                   // Rmin
                                          5.,                   // Rmax
                                          35.);                 // half length
    TGeoCombiTrans* t_part2 = new TGeoCombiTrans("t_part2", 0., 0., 65, fRefRot);
    t_part2->RegisterYourself();

    TGeoCompositeShape* pCBWorld = new TGeoCompositeShape(
        "Califa_box", " Califa_boxOut - (Califa_Centerpart1:t_part1  + Califa_Centerpart2:t_part2)");

    TGeoVolume* pWorld = new TGeoVolume("CalifaWorld", pCBWorld, pAirMedium);

    TGeoCombiTrans* t0 = new TGeoCombiTrans();
    TGeoCombiTrans* pGlobalc = GetGlobalPosition(t0);

    // add the sphere as Mother Volume
    top->AddNode(pWorld, 0, pGlobalc);

    // CALIFA CEPA USC VERSION (JUL 2020) PARAMETERS
    const Int_t N_ALV_TYPES = 3; // alveolar structures
    const Int_t N_CRY_TYPES = 12; // crystal elements
    const Double_t crystal_length = 20.0; //20cm

    ifstream in1, in2;
    in1.open("./CLF-CEPA_USC-onePart.txt");
    //in2.open("./CLF-CEPA_USC-oneCrystal.txt");
    Int_t counter = 0;
    Float_t x, y, z;

    Double_t wrapping_thickness = 0.0065; // in cm.
    Double_t cf_thickness = 0.0150; // in cm.
    // target reference in mm. OFFSET FROM UVIGO
    //TVector3 target_ref(4.1, 2304.0809, 325.0);
    TVector3 target_ref(0.0, 0.0, 0.0);

    // 3 geometries, 8 vertices, outer and inner: (3*8*2)
    TVector3 file_points[8];
    TVector3 points[N_ALV_TYPES * 8];
    TVector3 points_local[N_ALV_TYPES * 8];
    TVector3 points_inn[N_ALV_TYPES * 8];
    TVector3 points_inn_local[N_ALV_TYPES * 8];

    // 12 geometries, 8 vertices
    TVector3 points_cry[N_CRY_TYPES * 8];
    TVector3 points_cry_local[N_CRY_TYPES * 8];

    while (1)
    { // reading the file with all alveoli vertices
        in1 >> x >> y >> z;
        if (!in1.good())
            break;
        file_points[counter].SetXYZ(
            (x - target_ref.X()) / 10, (y - target_ref.Y()) / 10, (z - target_ref.Z()) / 10); // in cm;
        if (file_points[counter].X() > 100 || file_points[counter].Y() > 100 || file_points[counter].Z() > 100)
            cout << "WARNING: points exceed top volume!!" << endl;
        counter++;
        //printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    }
    if (counter != 8)
        cout << "PROBLEM! Counter=" << counter << endl;
    counter = 0;

    // external alveoli corners (calculated from external CEPA measurement)
    for(Int_t uOrD=0; uOrD<2; uOrD++)
    { //upper or lower face
        points[uOrD*4+0] = file_points[uOrD*4+0];       //crystal 1
        points[uOrD*4+1] = 0.5*(file_points[uOrD*4+0]+file_points[uOrD*4+1]);
        points[uOrD*4+2] = 0.25*(file_points[uOrD*4+0]+file_points[uOrD*4+1]+file_points[uOrD*4+2]+file_points[uOrD*4+3]);
        points[uOrD*4+3] = 0.5*(file_points[uOrD*4+3]+file_points[uOrD*4+0]);

        points[uOrD*4+8] = 0.5*(file_points[uOrD*4+0]+file_points[uOrD*4+1]);  //crystal 2
        points[uOrD*4+9] = file_points[uOrD*4+1];
        points[uOrD*4+10] = file_points[uOrD*4+2];
        points[uOrD*4+11] = 0.5*(file_points[uOrD*4+2]+file_points[uOrD*4+3]);

        points[uOrD*4+16] = 0.5*(file_points[uOrD*4+3]+file_points[uOrD*4+0]);  //crystal 4
        points[uOrD*4+17] = 0.25*(file_points[uOrD*4+0]+file_points[uOrD*4+1]+file_points[uOrD*4+2]+file_points[uOrD*4+3]);
        points[uOrD*4+18] = 0.5*(file_points[uOrD*4+2]+file_points[uOrD*4+3]);
        points[uOrD*4+19] = file_points[uOrD*4+3];
    }

    // centers of faces
    TVector3 center[N_ALV_TYPES * 2];     // 3 geometries, 2 face centers, (3*2)
    TVector3 x_uni[N_ALV_TYPES * 2];      // unit vectors for each face
    TVector3 y_uni[N_ALV_TYPES * 2];
    TVector3 z_uni[N_ALV_TYPES * 2];
    TVector3 center_inn[N_ALV_TYPES * 2];     // 3 geometries, 2 face centers, (3*2)
    TVector3 x_inn_uni[N_ALV_TYPES * 2];      // unit vectors for each face
    TVector3 y_inn_uni[N_ALV_TYPES * 2];
    TVector3 z_inn_uni[N_ALV_TYPES * 2];
    TVector3 center_cry[N_CRY_TYPES * 2]; // 12 types of crystals, 2 face centers (12*2)
    TVector3 x_uni_cry[N_CRY_TYPES * 2];  // unit vectors for each face for crystals
    TVector3 y_uni_cry[N_CRY_TYPES * 2];
    TVector3 z_uni_cry[N_CRY_TYPES * 2];
    TRotation rot[N_ALV_TYPES * 2];       // calculated in each face, but only 3 are really different if all is ok
    TRotation rot_inn[N_ALV_TYPES * 2];   // calculated in each face, but only 3 are really different if all is ok
    TRotation rot_cry[N_CRY_TYPES * 2];   // only a few are really different if all is ok
    // volume centers
    TVector3 alv_cm[N_ALV_TYPES];         // 3 geometries
    TVector3 alv_cm_rot[N_ALV_TYPES];     // 3 geometries, after final rotation
    TVector3 alv_inn_cm[N_ALV_TYPES];     // 3 geometries
    TVector3 alv_inn_cm_rot[N_ALV_TYPES]; // 3 geometries, after final rotation
    TVector3 cry_cm[N_CRY_TYPES];         // 12 types of crystals

    // The center of the faces are first calculated. Then, the unit vectors defining the axis in each faces
    // Third, the rotation moving from the lab system to the unit vectors previously found. To define the
    // volume in Arb8 style, we need the 8 corners in the local frustrum coordinates. Then, we should express
    // the vertices in the coordinate system of the volume center of mass (cm)
    for (Int_t i = 0; i < N_ALV_TYPES * 2; i++)
    { // for 3 geometries, 2 face centers (3*2)
        center[i] = points[i * 4] + points[i * 4 + 1] + points[i * 4 + 2] + points[i * 4 + 3]; // face centers
        center[i] *= 0.25;                                                                     // face centers
        // center[i].Print();
        z_uni[i] = (points[i * 4 + 1] - points[i * 4]).Cross(points[i * 4 + 2] - points[i * 4 + 1]);
        z_uni[i] = z_uni[i].Unit(); // normal to face center
        x_uni[i] = points[i * 4 + 2] - points[i * 4 + 1]; //MODIFIED FROM BARREL+IPHOS DEFINITION!!!!

        x_uni[i] = x_uni[i].Unit();          // unit along X
        y_uni[i] = z_uni[i].Cross(x_uni[i]); // unit along Y
        // x_uni[i].Print();  y_uni[i].Print();   z_uni[i].Print();
        // calculate rotation matrix for the 3 geometries (should be repeated 4 times, just checking)
        rot[i].SetZAxis(z_uni[i], x_uni[i]);
    }
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    { // for 3 geometries
        alv_cm[i] = center[2 * i] + center[2 * i + 1];
        alv_cm[i] *= 0.5; // volume center for all alv (outer and inner)
        // alv_cm[i].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each alveolus
            points_local[i * 8 + j] = rot[2 * i].Inverse() * (points[i * 8 + j] - alv_cm[i]);
            //cout<< "Points in cm coordinates: "<< endl;
            //points_local[i*8+j].Print();
        }
    }

/*
    cout<< "Points local in cm coordinates: "<< endl;
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
        for (Int_t j = 0; j < 8; j++)
             points_local[i*8+j].Print();
*/

    // same for inner points
    TVector3 d0(cf_thickness,-cf_thickness,cf_thickness);
    TVector3 d1(cf_thickness,cf_thickness,cf_thickness);
    TVector3 d2(-cf_thickness,cf_thickness,cf_thickness);
    TVector3 d3(-cf_thickness,-cf_thickness,cf_thickness);
    TVector3 d4(cf_thickness,-cf_thickness,-cf_thickness);
    TVector3 d5(cf_thickness,cf_thickness,-cf_thickness);
    TVector3 d6(-cf_thickness,cf_thickness,-cf_thickness);
    TVector3 d7(-cf_thickness,-cf_thickness,-cf_thickness);
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    { // for 3 geometries
      points_inn_local[i * 8 + 0] = points_local[i * 8 + 0] + d0;
      points_inn_local[i * 8 + 1] = points_local[i * 8 + 1] + d1;
      points_inn_local[i * 8 + 2] = points_local[i * 8 + 2] + d2;
      points_inn_local[i * 8 + 3] = points_local[i * 8 + 3] + d3;
      points_inn_local[i * 8 + 4] = points_local[i * 8 + 4] + d4;
      points_inn_local[i * 8 + 5] = points_local[i * 8 + 5] + d5;
      points_inn_local[i * 8 + 6] = points_local[i * 8 + 6] + d6;
      points_inn_local[i * 8 + 7] = points_local[i * 8 + 7] + d7;
      //for (Int_t j = 0; j < 8; j++) points_inn_local[i*8+j].Print();
      for (Int_t j = 0; j < 8; j++)
      { // for the 8 vertices of each inner alveolus (inverse of the normal point to point local)
          points_inn[i * 8 + j] = alv_cm[i] + rot[2 * i] * points_inn_local[i * 8 + j];
          //cout<< "Points in cm coordinates: "<< endl;
          //points_inn[i*8+j].Print();
      }
    }

    for (Int_t i = 0; i < N_ALV_TYPES * 2; i++)
    { // for 3 geometries, 2 face centers (3*2)
        center_inn[i] = points_inn[i * 4] + points_inn[i * 4 + 1] + points_inn[i * 4 + 2] + points_inn[i * 4 + 3]; // face centers
        center_inn[i] *= 0.25;                                                                     // face centers
        // center[i].Print();
        z_inn_uni[i] = (points_inn[i * 4 + 1] - points_inn[i * 4]).Cross(points_inn[i * 4 + 2] - points_inn[i * 4 + 1]);
        z_inn_uni[i] = z_inn_uni[i].Unit(); // normal to face center
        x_inn_uni[i] = points_inn[i * 4 + 2] - points_inn[i * 4 + 1]; //MODIFIED FROM BARREL+IPHOS DEFINITION!!!!

        x_inn_uni[i] = x_inn_uni[i].Unit();              // unit along X
        y_inn_uni[i] = z_inn_uni[i].Cross(x_inn_uni[i]); // unit along Y
        // x_uni[i].Print();  y_uni[i].Print();   z_uni[i].Print();
        // calculate rotation matrix for the 3 geometries (should be repeated 4 times, just checking)
        rot_inn[i].SetZAxis(z_inn_uni[i], x_inn_uni[i]);
    }
    // inner points_cry
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    { // for 3 geometries
        alv_inn_cm[i] = center_inn[2 * i] + center_inn[2 * i + 1];
        alv_inn_cm[i] *= 0.5; // volume center for all alv (outer and inner)
        // alv_cm[i].Print();
        //NEXT LINES ARE ONLY FOR CHECKING THE ROTATION
        //for (Int_t j = 0; j < 8; j++)
        //{ // for the 8 vertices of each alveolus
            //points_inn_local[i * 8 + j] = rot_inn[2 * i].Inverse() * (points_inn[i * 8 + j] - alv_inn_cm[i]);
            // cout<< "Points in cm coordinates: "<< endl;
            // points_local[i*8+j].Print();
        //}
    }

    // crystals corners
    Int_t val=0;
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    { //four crystals per alveoli, points calculated for each face
      for(Int_t uOrD=0; uOrD<2; uOrD++){ //upper or lower face
          val=i*8+uOrD*4; //the running index
          points_cry[i*32+uOrD*4+0] = points_inn[val+0];       //crystal 1
          points_cry[i*32+uOrD*4+1] = 0.5*(points_inn[val+0]+points_inn[val+1]);
          points_cry[i*32+uOrD*4+2] = 0.25*(points_inn[val+0]+points_inn[val+1]+points_inn[val+2]+points_inn[val+3]);
          points_cry[i*32+uOrD*4+3] = 0.5*(points_inn[val+3]+points_inn[val+0]);

          points_cry[i*32+uOrD*4+8] = 0.5*(points_inn[val+0]+points_inn[val+1]);  //crystal 2
          points_cry[i*32+uOrD*4+9] = points_inn[val+1];
          points_cry[i*32+uOrD*4+10] = 0.5*(points_inn[val+1]+points_inn[val+2]);
          points_cry[i*32+uOrD*4+11] = 0.25*(points_inn[val+0]+points_inn[val+1]+points_inn[val+2]+points_inn[val+3]);

          points_cry[i*32+uOrD*4+16] = 0.25*(points_inn[val+0]+points_inn[val+1]+points_inn[val+2]+points_inn[val+3]); //crystal 3
          points_cry[i*32+uOrD*4+17] = 0.5*(points_inn[val+1]+points_inn[val+2]);
          points_cry[i*32+uOrD*4+18] = points_inn[val+2];
          points_cry[i*32+uOrD*4+19] = 0.5*(points_inn[val+2]+points_inn[val+3]);

          points_cry[i*32+uOrD*4+24] = 0.5*(points_inn[val+3]+points_inn[val+0]);  //crystal 4
          points_cry[i*32+uOrD*4+25] = 0.25*(points_inn[val+0]+points_inn[val+1]+points_inn[val+2]+points_inn[val+3]);
          points_cry[i*32+uOrD*4+26] = 0.5*(points_inn[val+2]+points_inn[val+3]);
          points_cry[i*32+uOrD*4+27] = points_inn[val+3];
      }
    }

    // same procedure for crystals
    for (Int_t i = 0; i < N_CRY_TYPES * 2; i++)
    { // 12 types of crystals, 2 face centers (12*2)
        center_cry[i] =
            points_cry[i * 4] + points_cry[i * 4 + 1] + points_cry[i * 4 + 2] + points_cry[i * 4 + 3]; // face centers
        center_cry[i] *= 0.25;                                                                         // face centers
        // center_cry[i].Print();
        z_uni_cry[i] = (points_cry[i * 4 + 1] - points_cry[i * 4]).Cross(points_cry[i * 4 + 2] - points_cry[i * 4 + 1]);
        z_uni_cry[i] = z_uni_cry[i].Unit(); // normal to face center
        x_uni_cry[i] = points_cry[i * 4 + 2] - points_cry[i * 4 + 1]; //MODIFIED FROM BARREL+IPHOS DEFINITION!!!!

        x_uni_cry[i] = x_uni_cry[i].Unit();              // unit along X
        y_uni_cry[i] = z_uni_cry[i].Cross(x_uni_cry[i]); // unit along Y
        // x_uni_cry[i].Print();  y_uni_cry[i].Print();    z_uni_cry[i].Print();

        // calculate rotation matrix for the 85 types of crystals (should be repeated 4 times, just checking)
        rot_cry[i].SetZAxis(z_uni_cry[i], x_uni_cry[i]);
    }
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    { // 12 types of crystals
        cry_cm[i] = center_cry[2 * i] + center_cry[2 * i + 1];
        cry_cm[i] *= 0.5; // volume center for each cry
        // cry_cm[i].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each crystal
            points_cry_local[i * 8 + j] = rot_cry[2 * i].Inverse() * (points_cry[i * 8 + j] - cry_cm[i]);
        }
    }

    //Make crystals of a given length... calculation in local system coordinates
    Double_t proj_crystal_length = 0;
    TVector3 side_vec[N_CRY_TYPES * 4];
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    {
      for (Int_t j = 0; j < 4; j++)
      {
        side_vec[i*4+j] = points_cry_local[i*8+j+4]-points_cry_local[i*8+j];
        proj_crystal_length = crystal_length / (side_vec[i*4+j].Unit()).Z();
        points_cry_local[i*8+j+4] = points_cry_local[i*8+j] + proj_crystal_length * side_vec[i*4+j].Unit();
      }
    }
    //Recalculating crystal vertex after length scaling
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
      for (Int_t j = 0; j < 8; j++)
        points_cry[i * 8 + j] = cry_cm[i] + rot_cry[2 * i] * points_cry_local[i * 8 + j];

    // after the modification of the crystal lenght, centers should be recalculated
    for (Int_t i = 0; i < N_CRY_TYPES * 2; i++)
    { // 12 types of crystals, 2 face centers (12*2)
        center_cry[i] =
            points_cry[i * 4] + points_cry[i * 4 + 1] + points_cry[i * 4 + 2] + points_cry[i * 4 + 3]; // face centers
        center_cry[i] *= 0.25;                                                                         // face centers
        // center_cry[i].Print();
        z_uni_cry[i] = (points_cry[i * 4 + 1] - points_cry[i * 4]).Cross(points_cry[i * 4 + 2] - points_cry[i * 4 + 1]);
        z_uni_cry[i] = z_uni_cry[i].Unit(); // normal to face center
        x_uni_cry[i] = points_cry[i * 4 + 2] - points_cry[i * 4 + 1]; //MODIFIED FROM BARREL+IPHOS DEFINITION!!!!

        x_uni_cry[i] = x_uni_cry[i].Unit();              // unit along X
        y_uni_cry[i] = z_uni_cry[i].Cross(x_uni_cry[i]); // unit along Y
        // x_uni_cry[i].Print();  y_uni_cry[i].Print();    z_uni_cry[i].Print();

        // calculate rotation matrix for the 85 types of crystals (should be repeated 4 times, just checking)
        rot_cry[i].SetZAxis(z_uni_cry[i], x_uni_cry[i]);
    }
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    { // 12 types of crystals
        cry_cm[i] = center_cry[2 * i] + center_cry[2 * i + 1];
        cry_cm[i] *= 0.5; // volume center for each cry
        // cry_cm[i].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each crystal
            points_cry_local[i * 8 + j] = rot_cry[2 * i].Inverse() * (points_cry[i * 8 + j] - cry_cm[i]);
        }
    }

    // location of the crystals in the alveoli
    TVector3 cry_position[N_CRY_TYPES];       // 12 types of crystals
    TVector3 cry_position_local[N_CRY_TYPES]; // 12 types of crystals

    // Relative Crystal rotation in each alveoli. Obtained from the crystal unit vector in
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    {// four crystals per alv
      cry_position[i] = cry_cm[i] - alv_cm[(Int_t)(i / 4)];
      cry_position_local[i] =
        rot[2*(Int_t)(i / 4)].Inverse() * (cry_cm[i] - alv_cm[(Int_t)(i / 4)]);
    }

        /*
        cout<< "Points local in cm coordinates: "<< endl;
        for (Int_t i = 0; i < N_ALV_TYPES; i++)
            for (Int_t j = 0; j < 8; j++)
                 points_local[i*8+j].Print();
        cout<< "Points inn local in cm coordinates: "<< endl;
        for (Int_t i = 0; i < N_ALV_TYPES; i++)
            for (Int_t j = 0; j < 8; j++)
                points_inn_local[i*8+j].Print();
        cout<< "Points (check) cry local in cm coordinates after reduction: "<< endl;
        for (Int_t i = 0; i < N_CRY_TYPES; i++)
            for (Int_t j = 0; j < 8; j++)
                points_cry_local[i * 8 + j].Print();
        */


    // Redefinition of vertices for the construction of the Alveoli, using TGeoArb8
    Double_t* vertices_Alv[N_ALV_TYPES]; // 3 geometries
    Double_t* vertices_inner_Alv[N_ALV_TYPES];
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        vertices_Alv[i] = new Double_t[16];
        vertices_inner_Alv[i] = new Double_t[16];
    }
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        for (Int_t j = 0; j < 8; j++)
        {  // reversing order for being clockwise filling TGeoArb8
          if(j==1 || j== 5){
            vertices_Alv[i][2 * j] = points_local[8 * i + j + 2].X();
            vertices_Alv[i][2 * j + 1] = points_local[8 * i + j + 2].Y();
            vertices_inner_Alv[i][2 * j] = points_inn_local[8 * i + j + 2].X();
            vertices_inner_Alv[i][2 * j + 1] = points_inn_local[8 * i + j + 2].Y();

          }
          else if(j==3 || j== 7){
            vertices_Alv[i][2 * j] = points_local[8 * i + j - 2].X();
            vertices_Alv[i][2 * j + 1] = points_local[8 * i + j - 2].Y();
            vertices_inner_Alv[i][2 * j] = points_inn_local[8 * i + j - 2].X();
            vertices_inner_Alv[i][2 * j + 1] = points_inn_local[8 * i + j - 2].Y();

          }
          else {
              vertices_Alv[i][2 * j] = points_local[8 * i + j].X();
              vertices_Alv[i][2 * j + 1] = points_local[8 * i + j].Y();
              vertices_inner_Alv[i][2 * j] = points_inn_local[8 * i + j].X();
              vertices_inner_Alv[i][2 * j + 1] = points_inn_local[8 * i + j].Y();
            }
        }
    }

/*
    cout<< "Points local in cm coordinates reordered: "<< endl;
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        for (Int_t j = 0; j < 16; j++)
        {
        cout << vertices_Alv[i][j] << endl;
        }
    }
*/

    // Redefinition of vertices for the construction of the Crystals, using TGeoArb8
    Double_t* vertices_Cry[N_CRY_TYPES];      // 96/8=12 (12 crystal types)
    Double_t* vertices_Cry_Wrap[N_CRY_TYPES]; // 96/8=12 (12 crystal types)
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    {
        vertices_Cry[i] = new Double_t[16];
        vertices_Cry_Wrap[i] = new Double_t[16];
    }
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    {
        for (Int_t j = 0; j < 8; j++)
        {
            if ((4 - j) > 0)
            { // reversing order for being clockwise filling TGeoArb8
                vertices_Cry_Wrap[i][2 * j] = points_cry_local[8 * i + 3 - j].X();
                vertices_Cry_Wrap[i][2 * j + 1] = points_cry_local[8 * i + 3 - j].Y();
            }
            else
            {
                vertices_Cry_Wrap[i][2 * j] = points_cry_local[8 * i + 11 - j].X();
                vertices_Cry_Wrap[i][2 * j + 1] = points_cry_local[8 * i + 11 - j].Y();
            }
        }
    }
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    {
        for (Int_t j = 0; j < 16; j++)
        {
            if (vertices_Cry_Wrap[i][j] > 0)
                vertices_Cry[i][j] = vertices_Cry_Wrap[i][j] - wrapping_thickness;
            else
                vertices_Cry[i][j] = vertices_Cry_Wrap[i][j] + wrapping_thickness;
        }
    }

    TGeoVolume** Alv_vol;
    Alv_vol = new TGeoVolume*[N_ALV_TYPES];
    TGeoVolume** Alv_inner_vol;
    Alv_inner_vol = new TGeoVolume*[N_ALV_TYPES];
    TGeoVolume** Cry_vol_wrap;
    Cry_vol_wrap = new TGeoVolume*[N_CRY_TYPES];
    TGeoVolume** Cry_vol;
    Cry_vol = new TGeoVolume*[N_CRY_TYPES];

    TString AlvGlobalName = "Alveolus_CUSC_";
    TString AlvGlobalNameInner = "InnerAlv_CUSC_";
    // Substitute names in previous array (CAD names) to simplify the R3BRoot code
    TString name_Alv[N_ALV_TYPES] = { "01", "02", "03" };
    Double_t halfLengthAlv[N_ALV_TYPES] = { points_local[4].Z(),
                                            points_local[12].Z(),
                                            points_local[20].Z() }; // cm
    Double_t halfLengthAlv_inner[N_ALV_TYPES] = { points_inn_local[4].Z(),
                                                  points_inn_local[12].Z(),
                                                  points_inn_local[20].Z() }; // cm
    TString WrapCryGlobalName = "WrapCry_CUSC_";
    TString CryGlobalName = "Crystal_CUSC_";
    TString name_Cry[4] = { "_1", "_2", "_3", "_4" };
    // For the moment same length as the inner alveoli
    Double_t halfLengthCry[N_ALV_TYPES] = { crystal_length/2,
                                            crystal_length/2,
                                            crystal_length/2 }; // cm

    TGeoRotation* rotUni = new TGeoRotation();
    TGeoRotation** rotAlv = new TGeoRotation*[N_ALV_TYPES];
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
        rotAlv[i] = new TGeoRotation();
    TGeoRotation** rotCry = new TGeoRotation*[N_CRY_TYPES];
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
        rotCry[i] = new TGeoRotation();
    Double_t rotEle[9];

    // rotation
    TGeoRotation** rotOnZ = new TGeoRotation*[8];
    for (Int_t i = 0; i < 8; i++)
    {
        rotOnZ[i] = new TGeoRotation();
    }
    for (Int_t i = 0; i < 8; i++)
    {
        rotOnZ[i]->RotateZ(-45.0 * i + 67.5); //67.5 is the offset to put the first alveoli below the first of the barrel
        //rotOnZ[i]->Print();
    }
    TRotation** rotationOnZ = new TRotation*[8];
    for (Int_t i = 0; i < 8; i++)
    {
        rotationOnZ[i] = new TRotation();
    }
    for (Int_t i = 0; i < 8; i++)
    {
        rotationOnZ[i]->RotateZ((i * -45.0 + 67.5) * TMath::Pi() / 180);
    }

    TGeoRotation* oR = new TGeoRotation();
    //oR->Print();
    oR->RotateY(-90);
    //oR->Print();
    TRotation* oR2 = new TRotation();
    //oR->Print();
    oR2->RotateY(-90* TMath::Pi() / 180);
    //cout << oR2->XX() << " " << oR2->XY() << " " << oR2->XZ() << endl;
    //cout << oR2->YX() << " " << oR2->YY() << " " << oR2->YZ() << endl;
    //cout << oR2->ZX() << " " << oR2->ZY() << " " << oR2->ZZ() << endl;

    TGeoRotation** rotAlvFinal = new TGeoRotation*[8 * N_ALV_TYPES];
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        for (Int_t j = 0; j < 8; j++)
        {
          //rotAlvFinal[i * 8 + j] = new TGeoRotation( (*rotOnZ[j]) * (*rotAlv[i]));
          rotAlvFinal[i * 8 + j] = new TGeoRotation( (*rotOnZ[j]) *  (*rotAlv[i]) * (*oR));
            //if(i==0) rotAlvFinal[i * 8 + j]->Print();
        }
    }

    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        Alv_vol[i] =
            gGeoManager->MakeArb8(AlvGlobalName + name_Alv[i], pCarbonFibreMedium, halfLengthAlv[i], vertices_Alv[i]);
        Alv_vol[i]->SetLineColor(kBlue);
        Alv_vol[i]->SetVisLeaves(kTRUE);
        Alv_vol[i]->SetVisibility(kTRUE);
        Alv_vol[i]->SetVisContainers(kTRUE);

        // TODO: CHANGE VACUUM TO AIR!!!
        Alv_inner_vol[i] = gGeoManager->MakeArb8(
            AlvGlobalNameInner + name_Alv[i], pAirMedium, halfLengthAlv_inner[i], vertices_inner_Alv[i]);
        Alv_inner_vol[i]->SetLineColor(kRed);
        Alv_inner_vol[i]->SetVisLeaves(kTRUE);
        Alv_inner_vol[i]->SetVisibility(kTRUE);
        Alv_inner_vol[i]->SetVisContainers(kTRUE);

        // four crystals per alv
        for (Int_t j = 0; j < 4; j++)
        {
            Cry_vol[4 * i + j] = gGeoManager->MakeArb8(CryGlobalName + name_Alv[i] + name_Cry[j],
                                                           pCsIMedium,
                                                           halfLengthCry[i] - wrapping_thickness,
                                                           vertices_Cry[4 * i + j]);
            Cry_vol[4 * i + j]->SetLineColor(kMagenta);
            Cry_vol[4 * i + j]->SetVisLeaves(kTRUE);
            Cry_vol[4 * i + j]->SetVisibility(kTRUE);
            Cry_vol[4 * i + j]->SetVisContainers(kTRUE);

            Cry_vol_wrap[4 * i + j] = gGeoManager->MakeArb8(WrapCryGlobalName + name_Alv[i] + name_Cry[j],
                                                                pWrappingMedium,
                                                                halfLengthCry[i],
                                                                vertices_Cry_Wrap[4 * i + j]);
            Cry_vol_wrap[4 * i + j]->SetLineColor(kGreen);
            Cry_vol_wrap[4 * i + j]->SetVisLeaves(kTRUE);
            Cry_vol_wrap[4 * i + j]->SetVisibility(kTRUE);
            Cry_vol_wrap[4 * i + j]->SetVisContainers(kTRUE);

            Cry_vol_wrap[4 * i + j]->AddNode(Cry_vol[4 * i + j], 0, new TGeoCombiTrans(0, 0, 0, rotUni));
            Alv_inner_vol[i]->AddNode(Cry_vol_wrap[4 * i + j],
                                      0,
                                      new TGeoCombiTrans(cry_position_local[4 * i + j].X(),
                                                         cry_position_local[4 * i + j].Y(),
                                                         cry_position_local[4 * i + j].Z(),
                                                         rotUni));
        }

        // Inner volume center is displaced 150 microns along Z
        Alv_vol[i]->AddNode(Alv_inner_vol[i], 0, new TGeoCombiTrans(0, 0, cf_thickness, rotUni));

        rotEle[0] = rot[2 * i].XX();
        rotEle[1] = rot[2 * i].XY();
        rotEle[2] = rot[2 * i].XZ();
        rotEle[3] = rot[2 * i].YX();
        rotEle[4] = rot[2 * i].YY();
        rotEle[5] = rot[2 * i].YZ();
        rotEle[6] = rot[2 * i].ZX();
        rotEle[7] = rot[2 * i].ZY();
        rotEle[8] = rot[2 * i].ZZ();
        rotAlv[i]->SetMatrix(rotEle);

        for (Int_t j = 0; j < 8; j++)
        { // rotation around Z
            // rotAlvFinal[i * 32 + j] = new TGeoRotation((*rotOnZ[j]) * (*rotAlv[i]));
            // alv_cm_rot[2 * i] = (*rotationOnZ[j]) * alv_cm[2 * i];

            rotAlvFinal[i * 8 + j] = new TGeoRotation((*rotOnZ[j]) * (*oR) * (*rotAlv[i]));
            alv_cm_rot[i] =  (*rotationOnZ[j]) * (*oR2) * alv_cm[i];
            //if (isCrystalInstalled(i + 1, j))
            //{
                pWorld->AddNode(Alv_vol[i],
                                j,
                                new TGeoCombiTrans(alv_cm_rot[i].X(),
                                                   alv_cm_rot[i].Y(),
                                                   alv_cm_rot[i].Z(),
                                                   rotAlvFinal[i * 8 + j])
                                                    );
            //}
        }
    }
    // gGeoManager->SetVisOption(0);
    // gGeoManager->SetVisLevel(4); //only vis the indicated level

    gGeoManager->CloseGeometry();
    // gGeoManager->CheckGeometryFull();
    gGeoManager->CheckOverlaps(0.0001);
    gGeoManager->PrintOverlaps();
    gGeoManager->Test();

    TFile* geoFile = new TFile(geoFileName, "RECREATE");
    top->Write();
    geoFile->Close();
    // --------------------------------------------------------------------------
}

TGeoCombiTrans* GetGlobalPosition(TGeoCombiTrans* fRef)
{
    if (fLocalTrans == kTRUE)
    {

        if ((fThetaX == 0) && (fThetaY == 0) && (fThetaZ == 0) && (fX == 0) && (fY == 0) && (fZ == 0))
            return fRef;

        // X axis
        Double_t xAxis[3] = { 1., 0., 0. };
        Double_t yAxis[3] = { 0., 1., 0. };
        Double_t zAxis[3] = { 0., 0., 1. };
        // Reference Rotation
        fRefRot = fRef->GetRotation();

        if (fRefRot)
        {
            Double_t mX[3] = { 0., 0., 0. };
            Double_t mY[3] = { 0., 0., 0. };
            Double_t mZ[3] = { 0., 0., 0. };

            fRefRot->LocalToMasterVect(xAxis, mX);
            fRefRot->LocalToMasterVect(yAxis, mY);
            fRefRot->LocalToMasterVect(zAxis, mZ);

            Double_t a[4] = { mX[0], mX[1], mX[2], fThetaX };
            Double_t b[4] = { mY[0], mY[1], mY[2], fThetaY };
            Double_t c[4] = { mZ[0], mZ[1], mZ[2], fThetaZ };

            ROOT::Math::AxisAngle aX(a, a + 4);
            ROOT::Math::AxisAngle aY(b, b + 4);
            ROOT::Math::AxisAngle aZ(c, c + 4);

            ROOT::Math::Rotation3D fMatX(aX);
            ROOT::Math::Rotation3D fMatY(aY);
            ROOT::Math::Rotation3D fMatZ(aZ);

            ROOT::Math::Rotation3D fRotXYZ = (fMatZ * (fMatY * fMatX));

            // cout << fRotXYZ << endl;

            Double_t fRotable[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            fRotXYZ.GetComponents(fRotable[0],
                                  fRotable[3],
                                  fRotable[6],
                                  fRotable[1],
                                  fRotable[4],
                                  fRotable[7],
                                  fRotable[2],
                                  fRotable[5],
                                  fRotable[8]);
            TGeoRotation* pRot = new TGeoRotation();
            pRot->SetMatrix(fRotable);
            TGeoCombiTrans* pTmp = new TGeoCombiTrans(*fGlobalTrans, *pRot);

            // ne peut pas etre applique ici
            // il faut differencier trans et rot dans la multi.
            TGeoRotation rot_id;
            rot_id.SetAngles(0.0, 0.0, 0.0);

            TGeoCombiTrans c1;
            c1.SetRotation(rot_id);
            const Double_t* t = pTmp->GetTranslation();
            c1.SetTranslation(t[0], t[1], t[2]);

            TGeoCombiTrans c2;
            c2.SetRotation(rot_id);
            const Double_t* tt = fRefRot->GetTranslation();
            c2.SetTranslation(tt[0], tt[1], tt[2]);

            TGeoCombiTrans cc = c1 * c2;

            TGeoCombiTrans c3;
            c3.SetRotation(pTmp->GetRotation());
            TGeoCombiTrans c4;
            c4.SetRotation(fRefRot);

            TGeoCombiTrans ccc = c3 * c4;

            TGeoCombiTrans pGlobal;
            pGlobal.SetRotation(ccc.GetRotation());
            const Double_t* allt = cc.GetTranslation();
            pGlobal.SetTranslation(allt[0], allt[1], allt[2]);

            return (new TGeoCombiTrans(pGlobal));
        }
        else
        {

            cout << "-E- R3BDetector::GetGlobalPosition() \
              	      No. Ref. Transformation defined ! "
                 << endl;
            cout << "-E- R3BDetector::GetGlobalPosition() \
              	      cannot create Local Transformation "
                 << endl;
            return NULL;
        } //! fRefRot
    }
    else
    {
        // Lab Transf.
        if ((fPhi == 0) && (fTheta == 0) && (fPsi == 0) && (fX == 0) && (fY == 0) && (fZ == 0))
            return fRef;

        return (new TGeoCombiTrans(*fGlobalTrans, *fGlobalRot));
    }
}

Bool_t isCrystalInstalled(Int_t alvType, Int_t alveolusCopy)
{
    //
    // reproduces partially the algorithm of R3BCalifaGeometry::GetCrystalId(const char* volumePath)
    //

    ifstream wc1;
    wc1.open("./califa_InstalledCrystals_Nov2019.txt");

    Bool_t found = kFALSE;
    /*
    Int_t crystalId = 0;
    Int_t counter = 0;
    Int_t read = 0;
    Int_t cryType = 1;             // the first crystal of the alveoli... if not present, alveoli is removed.
    Int_t installedCrystals[4864]; // 2* 2432 just to cope numbers for crystals with both amplifications
    for (Int_t i = 0; i < 4864; i++)
        installedCrystals[i] = 0;

    while (1)
    { // reading the file with all alveoli vertices
        wc1 >> read;
        if (!wc1.good())
            break;
        installedCrystals[counter] = read;
        counter++;
    }

    if (alvType == 1)
        crystalId = 1 + alveolusCopy; // first alveoli ring, one crystal per alveolus
    else if (alvType < 20)
        crystalId = 33 + (alvType - 2) * 128 + alveolusCopy * 4 + (cryType - 1); // four crystal per alveolus
    else
        crystalId = 2337 + (alvType - 20) * 24 + alveolusCopy * 3 + (cryType - 1); // three crystal per alveolus

    if (crystalId < 1 || crystalId > 2432)
    { // crystalId runs from 1 to 2432
        cout << "R3BCalifaGeometry: Wrong crystal number ";
        cout << "---- crystalId: " << crystalId << endl;
        return 0;
    }
    for (Int_t i = 0; i < 4864; i++)
    {
        if (crystalId == installedCrystals[i])
            found = kTRUE;
    }
    */
    return found;
}
