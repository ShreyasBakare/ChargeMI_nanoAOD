#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//#include <nano9Ana.h>
/*
This is a driver script.
It decides which code to run over which sample, the names
of output files and so on.
*/


void ana(int sample=0){
  const char *hstfilename, *sumfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events"); //"Events"
  //Declare an instance of our code class
  nano9Ana m_selec;
  
  if(sample==0){
    //Add one file to chain. This is the input file.
    TString dy = "inputs/hst_DY";
    TString root = ".root";
    TString name = dy + root;
    chain->Add(name);
    //Set Names of outputfiles
    hstfilename = "output/hst_DY.root";
    sumfilename = "output/sum_DY.txt";
    //Set some options
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2016);
  }
  if(sample==1){
    chain->Add("inputs/TTTo2L2Nu.root");
    hstfilename = "output/hst_tt.root";
    sumfilename = "output/sum_tt.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2016);
  }
  if(sample==2){
    //Add one file to chain. This is the input file.
    chain->Add("inputs/TTZ.root");
    //Set Names of outputfiles
    hstfilename = "output/hst_TTZ.root";
    sumfilename = "output/sum_TTZ.txt";
    //Set some options
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }
  if(sample==3){
    chain->Add("inputs/WZ.root");
    hstfilename = "output/hst_WZ.root";
    sumfilename = "output/sum_WZ.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }
  if(sample==4){
    chain->Add("inputs/ZZ.root");
    hstfilename = "output/hst_ZZ.root";
    sumfilename = "output/sum_ZZ.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2018);
  }
  
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  // Set some more options.. set the output file names.
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);

}

int main(int argc, char *argv[])
{
  if (argc<2){
    std::cout<<" please give one integer argument "<<std::endl;
    return 0;
  }
  int sample_id = atoi(argv[1]);

  ana(sample_id);
  return 0;

}
