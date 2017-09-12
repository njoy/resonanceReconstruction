#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "resonanceReconstruction.hpp"

using namespace njoy::resonanceReconstruction;
using namespace dimwits;

ENDFtk::section::Type<2>
resonances();

SCENARIO( "Integration test" ){
/*
  auto Rh105 = resonances();
  const auto& isotope = Rh105.isotopes().front();
  const auto& energyRange = isotope.energyRanges().front();
  const auto& slbw =
    std::get< ENDFtk::resolved::SingleLevelBreitWigner >( energyRange );
*/  
}

std::string Rhodium105Resonances();

ENDFtk::section::Type<2>
resonances(){
  auto string = Rhodium105Resonances.begin();
  auto start = string.begin()
  auto it = start;
  auto end = string.end();
  auto lineNumber = 1;

  auto head = ENDFtk::HEAD( it, end, lineNumber );
  return ENDFtk::section::Type<2>( head, it, end, lineNumber, 4531 );
}

std::string Rhodium105Resonances(){
  return
    " 4.510500+4 1.040000+2          0          0          1          04531 2151    1\n"
    " 4.510500+4 1.000000+0          0          0          2          04531 2151    2\n"
    " 1.000000-5 7.500000+0          1          1          0          04531 2151    3\n"
    " 5.000000-1 6.200000-1          0          0          1          04531 2151    4\n"
    " 1.040050+2 0.000000+0          0          0         12          24531 2151    5\n"
    "-5.000000+0 1.000000+0 1.610000+0 1.450000+0 1.600000-1 0.000000+04531 2151    6\n"
    " 5.000000+0 1.000000+0 4.900000-1 3.300000-1 1.600000-1 0.000000+04531 2151    7\n"
    " 7.500000+0 1.000000+5          2          2          0          04531 2151    8\n"
    " 3.500000+0 6.207500-1          0          0          3          04531 2151    9\n"
    " 1.040050+2 0.000000+0          0          0          2          04531 2151   10\n"
    " 3.000000+0 0.000000+0          5          0        174         284531 2151   11\n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04531 2151   12\n"
    " 7.500000+0 3.410700+1 0.000000+0 1.500700-3 1.500000-1 0.000000+04531 2151   13\n"
    " 1.000000+1 3.267800+1 0.000000+0 1.437800-3 1.500000-1 0.000000+04531 2151   14\n"
    " 1.000000+2 3.399100+1 0.000000+0 1.495600-3 1.500000-1 0.000000+04531 2151   15\n"
    " 2.000000+2 3.419800+1 0.000000+0 1.504700-3 1.500000-1 0.000000+04531 2151   16\n"
    " 3.000000+2 3.429900+1 0.000000+0 1.509200-3 1.500000-1 0.000000+04531 2151   17\n"
    " 4.000000+2 3.436700+1 0.000000+0 1.512100-3 1.500000-1 0.000000+04531 2151   18\n"
    " 5.000000+2 3.441100+1 0.000000+0 1.514100-3 1.500000-1 0.000000+04531 2151   19\n"
    " 6.000000+2 3.444800+1 0.000000+0 1.515700-3 1.500000-1 0.000000+04531 2151   20\n"
    " 8.000000+2 3.449900+1 0.000000+0 1.517900-3 1.500000-1 0.000000+04531 2151   21\n"
    " 1.000000+3 3.452900+1 0.000000+0 1.519300-3 1.500000-1 0.000000+04531 2151   22\n"
    " 1.500000+3 3.458800+1 0.000000+0 1.521900-3 1.500000-1 0.000000+04531 2151   23\n"
    " 2.000000+3 3.463000+1 0.000000+0 1.523700-3 1.500000-1 0.000000+04531 2151   24\n"
    " 3.000000+3 3.468400+1 0.000000+0 1.526100-3 1.500000-1 0.000000+04531 2151   25\n"
    " 4.000000+3 3.457000+1 0.000000+0 1.521100-3 1.500000-1 0.000000+04531 2151   26\n"
    " 5.000000+3 3.461500+1 0.000000+0 1.523100-3 1.500000-1 0.000000+04531 2151   27\n"
    " 6.000000+3 3.465700+1 0.000000+0 1.524900-3 1.500000-1 0.000000+04531 2151   28\n"
    " 7.000000+3 3.466600+1 0.000000+0 1.525300-3 1.500000-1 0.000000+04531 2151   29\n"
    " 8.000000+3 3.468700+1 0.000000+0 1.526200-3 1.500000-1 0.000000+04531 2151   30\n"
    " 1.000000+4 3.471400+1 0.000000+0 1.527400-3 1.500000-1 0.000000+04531 2151   31\n"
    " 1.500000+4 3.469900+1 0.000000+0 1.526800-3 1.500000-1 0.000000+04531 2151   32\n"
    " 2.500000+4 3.464300+1 0.000000+0 1.524300-3 1.500000-1 0.000000+04531 2151   33\n"
    " 3.000000+4 3.452100+1 0.000000+0 1.518900-3 1.500000-1 0.000000+04531 2151   34\n"
    " 4.000000+4 3.423700+1 0.000000+0 1.506400-3 1.500000-1 0.000000+04531 2151   35\n"
    " 5.000000+4 3.392600+1 0.000000+0 1.492800-3 1.500000-1 0.000000+04531 2151   36\n"
    " 6.000000+4 3.359200+1 0.000000+0 1.478000-3 1.500000-1 0.000000+04531 2151   37\n"
    " 7.000000+4 3.326400+1 0.000000+0 1.463600-3 1.500000-1 0.000000+04531 2151   38\n"
    " 8.000000+4 3.294100+1 0.000000+0 1.449400-3 1.500000-1 0.000000+04531 2151   39\n"
    " 1.000000+5 3.231300+1 0.000000+0 1.421800-3 1.500000-1 0.000000+04531 2151   40\n"
    " 4.000000+0 0.000000+0          5          0        174         284531 2151   41\n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04531 2151   42\n"
    " 7.500000+0 2.652700+1 0.000000+0 1.167200-3 1.500000-1 0.000000+04531 2151   43\n"
    " 1.000000+1 2.541600+1 0.000000+0 1.118300-3 1.500000-1 0.000000+04531 2151   44\n"
    " 1.000000+2 2.643800+1 0.000000+0 1.163300-3 1.500000-1 0.000000+04531 2151   45\n"
    " 2.000000+2 2.659900+1 0.000000+0 1.170300-3 1.500000-1 0.000000+04531 2151   46\n"
    " 3.000000+2 2.667700+1 0.000000+0 1.173800-3 1.500000-1 0.000000+04531 2151   47\n"
    " 4.000000+2 2.673000+1 0.000000+0 1.176100-3 1.500000-1 0.000000+04531 2151   48\n"
    " 5.000000+2 2.676400+1 0.000000+0 1.177600-3 1.500000-1 0.000000+04531 2151   49\n"
    " 6.000000+2 2.679300+1 0.000000+0 1.178900-3 1.500000-1 0.000000+04531 2151   50\n"
    " 8.000000+2 2.683200+1 0.000000+0 1.180600-3 1.500000-1 0.000000+04531 2151   51\n"
    " 1.000000+3 2.685600+1 0.000000+0 1.181700-3 1.500000-1 0.000000+04531 2151   52\n"
    " 1.500000+3 2.690200+1 0.000000+0 1.183700-3 1.500000-1 0.000000+04531 2151   53\n"
    " 2.000000+3 2.693400+1 0.000000+0 1.185100-3 1.500000-1 0.000000+04531 2151   54\n"
    " 3.000000+3 2.697600+1 0.000000+0 1.187000-3 1.500000-1 0.000000+04531 2151   55\n"
    " 4.000000+3 2.688800+1 0.000000+0 1.183100-3 1.500000-1 0.000000+04531 2151   56\n"
    " 5.000000+3 2.692300+1 0.000000+0 1.184600-3 1.500000-1 0.000000+04531 2151   57\n"
    " 6.000000+3 2.695500+1 0.000000+0 1.186000-3 1.500000-1 0.000000+04531 2151   58\n"
    " 7.000000+3 2.696300+1 0.000000+0 1.186400-3 1.500000-1 0.000000+04531 2151   59\n"
    " 8.000000+3 2.697900+1 0.000000+0 1.187100-3 1.500000-1 0.000000+04531 2151   60\n"
    " 1.000000+4 2.700000+1 0.000000+0 1.188000-3 1.500000-1 0.000000+04531 2151   61\n"
    " 1.500000+4 2.698800+1 0.000000+0 1.187500-3 1.500000-1 0.000000+04531 2151   62\n"
    " 2.500000+4 2.694400+1 0.000000+0 1.185600-3 1.500000-1 0.000000+04531 2151   63\n"
    " 3.000000+4 2.684900+1 0.000000+0 1.181400-3 1.500000-1 0.000000+04531 2151   64\n"
    " 4.000000+4 2.662800+1 0.000000+0 1.171700-3 1.500000-1 0.000000+04531 2151   65\n"
    " 5.000000+4 2.638700+1 0.000000+0 1.161000-3 1.500000-1 0.000000+04531 2151   66\n"
    " 6.000000+4 2.612700+1 0.000000+0 1.149600-3 1.500000-1 0.000000+04531 2151   67\n"
    " 7.000000+4 2.587200+1 0.000000+0 1.138400-3 1.500000-1 0.000000+04531 2151   68\n"
    " 8.000000+4 2.562000+1 0.000000+0 1.127300-3 1.500000-1 0.000000+04531 2151   69\n"
    " 1.000000+5 2.513200+1 0.000000+0 1.105800-3 1.500000-1 0.000000+04531 2151   70\n"
    " 1.040050+2 0.000000+0          1          0          4          04531 2151   71\n"
    " 2.000000+0 0.000000+0          5          0        174         284531 2151   72\n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04531 2151   73\n"
    " 7.500000+0 4.774900+1 0.000000+0 1.957700-2 1.500000-1 0.000000+04531 2151   74\n"
    " 1.000000+1 4.574900+1 0.000000+0 1.875700-2 1.500000-1 0.000000+04531 2151   75\n"
    " 1.000000+2 4.758800+1 0.000000+0 1.951100-2 1.500000-1 0.000000+04531 2151   76\n"
    " 2.000000+2 4.787800+1 0.000000+0 1.963000-2 1.500000-1 0.000000+04531 2151   77\n"
    " 3.000000+2 4.801900+1 0.000000+0 1.968800-2 1.500000-1 0.000000+04531 2151   78\n"
    " 4.000000+2 4.811400+1 0.000000+0 1.972700-2 1.500000-1 0.000000+04531 2151   79\n"
    " 5.000000+2 4.817500+1 0.000000+0 1.975200-2 1.500000-1 0.000000+04531 2151   80\n"
    " 6.000000+2 4.822700+1 0.000000+0 1.977300-2 1.500000-1 0.000000+04531 2151   81\n"
    " 8.000000+2 4.829800+1 0.000000+0 1.980200-2 1.500000-1 0.000000+04531 2151   82\n"
    " 1.000000+3 4.834100+1 0.000000+0 1.982000-2 1.500000-1 0.000000+04531 2151   83\n"
    " 1.500000+3 4.842300+1 0.000000+0 1.985300-2 1.500000-1 0.000000+04531 2151   84\n"
    " 2.000000+3 4.848200+1 0.000000+0 1.987800-2 1.500000-1 0.000000+04531 2151   85\n"
    " 3.000000+3 4.855700+1 0.000000+0 1.990900-2 1.500000-1 0.000000+04531 2151   86\n"
    " 4.000000+3 4.839800+1 0.000000+0 1.984300-2 1.500000-1 0.000000+04531 2151   87\n"
    " 5.000000+3 4.846200+1 0.000000+0 1.986900-2 1.500000-1 0.000000+04531 2151   88\n"
    " 6.000000+3 4.851900+1 0.000000+0 1.989300-2 1.500000-1 0.000000+04531 2151   89\n"
    " 7.000000+3 4.853300+1 0.000000+0 1.989800-2 1.500000-1 0.000000+04531 2151   90\n"
    " 8.000000+3 4.856200+1 0.000000+0 1.991000-2 1.500000-1 0.000000+04531 2151   91\n"
    " 1.000000+4 4.860000+1 0.000000+0 1.992600-2 1.500000-1 0.000000+04531 2151   92\n"
    " 1.500000+4 4.857900+1 0.000000+0 1.991700-2 1.500000-1 0.000000+04531 2151   93\n"
    " 2.500000+4 4.850000+1 0.000000+0 1.988500-2 1.500000-1 0.000000+04531 2151   94\n"
    " 3.000000+4 4.832900+1 0.000000+0 1.981500-2 1.500000-1 0.000000+04531 2151   95\n"
    " 4.000000+4 4.793100+1 0.000000+0 1.965200-2 1.500000-1 0.000000+04531 2151   96\n"
    " 5.000000+4 4.749700+1 0.000000+0 1.947400-2 1.500000-1 0.000000+04531 2151   97\n"
    " 6.000000+4 4.702800+1 0.000000+0 1.928200-2 1.500000-1 0.000000+04531 2151   98\n"
    " 7.000000+4 4.657000+1 0.000000+0 1.909400-2 1.500000-1 0.000000+04531 2151   99\n"
    " 8.000000+4 4.611700+1 0.000000+0 1.890800-2 1.500000-1 0.000000+04531 2151  100\n"
    " 1.000000+5 4.523800+1 0.000000+0 1.854800-2 1.500000-1 0.000000+04531 2151  101\n"
    " 3.000000+0 0.000000+0          5          0        174         284531 2151  102\n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+04531 2151  103\n"
    " 7.500000+0 3.410700+1 0.000000+0 1.398400-2 1.500000-1 0.000000+04531 2151  104\n"
    " 1.000000+1 3.267800+1 0.000000+0 1.339800-2 1.500000-1 0.000000+04531 2151  105\n"
    " 1.000000+2 3.399100+1 0.000000+0 1.393600-2 1.500000-1 0.000000+04531 2151  106\n"
    " 2.000000+2 3.419800+1 0.000000+0 1.402100-2 1.500000-1 0.000000+04531 2151  107\n"
    " 3.000000+2 3.429900+1 0.000000+0 1.406300-2 1.500000-1 0.000000+04531 2151  108\n"
    " 4.000000+2 3.436700+1 0.000000+0 1.409000-2 1.500000-1 0.000000+04531 2151  109\n"
    " 5.000000+2 3.441100+1 0.000000+0 1.410800-2 1.500000-1 0.000000+04531 2151  110\n"
    " 6.000000+2 3.444800+1 0.000000+0 1.412400-2 1.500000-1 0.000000+04531 2151  111\n"
    " 8.000000+2 3.449900+1 0.000000+0 1.414400-2 1.500000-1 0.000000+04531 2151  112\n"
    " 1.000000+3 3.452900+1 0.000000+0 1.415700-2 1.500000-1 0.000000+04531 2151  113\n"
    " 1.500000+3 3.458800+1 0.000000+0 1.418100-2 1.500000-1 0.000000+04531 2151  114\n"
    " 2.000000+3 3.463000+1 0.000000+0 1.419800-2 1.500000-1 0.000000+04531 2151  115\n"
    " 3.000000+3 3.468400+1 0.000000+0 1.422000-2 1.500000-1 0.000000+04531 2151  116\n"
    " 4.000000+3 3.457000+1 0.000000+0 1.417400-2 1.500000-1 0.000000+04531 2151  117\n"
    " 5.000000+3 3.461500+1 0.000000+0 1.419200-2 1.500000-1 0.000000+04531 2151  118\n"
    " 6.000000+3 3.465700+1 0.000000+0 1.420900-2 1.500000-1 0.000000+04531 2151  119\n"
    " 7.000000+3 3.466600+1 0.000000+0 1.421300-2 1.500000-1 0.000000+04531 2151  120\n"
    " 8.000000+3 3.468700+1 0.000000+0 1.422200-2 1.500000-1 0.000000+04531 2151  121\n"
    " 1.000000+4 3.471400+1 0.000000+0 1.423300-2 1.500000-1 0.000000+04531 2151  122\n"
    " 1.500000+4 3.469900+1 0.000000+0 1.422700-2 1.500000-1 0.000000+04531 2151  123\n"
    " 2.500000+4 3.464300+1 0.000000+0 1.420400-2 1.500000-1 0.000000+04531 2151  124\n"
    " 3.000000+4 3.452100+1 0.000000+0 1.415300-2 1.500000-1 0.000000+04531 2151  125\n"
    " 4.000000+4 3.423700+1 0.000000+0 1.403700-2 1.500000-1 0.000000+04531 2151  126\n"
    " 5.000000+4 3.392600+1 0.000000+0 1.391000-2 1.500000-1 0.000000+04531 2151  127\n"
    " 6.000000+4 3.359200+1 0.000000+0 1.377300-2 1.500000-1 0.000000+04531 2151  128\n"
    " 7.000000+4 3.326400+1 0.000000+0 1.363800-2 1.500000-1 0.000000+04531 2151  129\n"
    " 8.000000+4 3.294100+1 0.000000+0 1.350600-2 1.500000-1 0.000000+04531 2151  130\n"
    " 1.000000+5 3.231300+1 0.000000+0 1.324800-2 1.500000-1 0.000000+04531 2151  131\n"
    " 4.000000+0 0.000000+0          5          0        174         284531 2151  132\n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+04531 2151  133\n"
    " 7.500000+0 2.652700+1 0.000000+0 1.087600-2 1.500000-1 0.000000+04531 2151  134\n"
    " 1.000000+1 2.541600+1 0.000000+0 1.042100-2 1.500000-1 0.000000+04531 2151  135\n"
    " 1.000000+2 2.643800+1 0.000000+0 1.083900-2 1.500000-1 0.000000+04531 2151  136\n"
    " 2.000000+2 2.659900+1 0.000000+0 1.090500-2 1.500000-1 0.000000+04531 2151  137\n"
    " 3.000000+2 2.667700+1 0.000000+0 1.093800-2 1.500000-1 0.000000+04531 2151  138\n"
    " 4.000000+2 2.673000+1 0.000000+0 1.095900-2 1.500000-1 0.000000+04531 2151  139\n"
    " 5.000000+2 2.676400+1 0.000000+0 1.097300-2 1.500000-1 0.000000+04531 2151  140\n"
    " 6.000000+2 2.679300+1 0.000000+0 1.098500-2 1.500000-1 0.000000+04531 2151  141\n"
    " 8.000000+2 2.683200+1 0.000000+0 1.100100-2 1.500000-1 0.000000+04531 2151  142\n"
    " 1.000000+3 2.685600+1 0.000000+0 1.101100-2 1.500000-1 0.000000+04531 2151  143\n"
    " 1.500000+3 2.690200+1 0.000000+0 1.103000-2 1.500000-1 0.000000+04531 2151  144\n"
    " 2.000000+3 2.693400+1 0.000000+0 1.104300-2 1.500000-1 0.000000+04531 2151  145\n"
    " 3.000000+3 2.697600+1 0.000000+0 1.106000-2 1.500000-1 0.000000+04531 2151  146\n"
    " 4.000000+3 2.688800+1 0.000000+0 1.102400-2 1.500000-1 0.000000+04531 2151  147\n"
    " 5.000000+3 2.692300+1 0.000000+0 1.103800-2 1.500000-1 0.000000+04531 2151  148\n"
    " 6.000000+3 2.695500+1 0.000000+0 1.105200-2 1.500000-1 0.000000+04531 2151  149\n"
    " 7.000000+3 2.696300+1 0.000000+0 1.105500-2 1.500000-1 0.000000+04531 2151  150\n"
    " 8.000000+3 2.697900+1 0.000000+0 1.106100-2 1.500000-1 0.000000+04531 2151  151\n"
    " 1.000000+4 2.700000+1 0.000000+0 1.107000-2 1.500000-1 0.000000+04531 2151  152\n"
    " 1.500000+4 2.698800+1 0.000000+0 1.106500-2 1.500000-1 0.000000+04531 2151  153\n"
    " 2.500000+4 2.694400+1 0.000000+0 1.104700-2 1.500000-1 0.000000+04531 2151  154\n"
    " 3.000000+4 2.684900+1 0.000000+0 1.100800-2 1.500000-1 0.000000+04531 2151  155\n"
    " 4.000000+4 2.662800+1 0.000000+0 1.091800-2 1.500000-1 0.000000+04531 2151  156\n"
    " 5.000000+4 2.638700+1 0.000000+0 1.081900-2 1.500000-1 0.000000+04531 2151  157\n"
    " 6.000000+4 2.612700+1 0.000000+0 1.071200-2 1.500000-1 0.000000+04531 2151  158\n"
    " 7.000000+4 2.587200+1 0.000000+0 1.060800-2 1.500000-1 0.000000+04531 2151  159\n"
    " 8.000000+4 2.562000+1 0.000000+0 1.050400-2 1.500000-1 0.000000+04531 2151  160\n"
    " 1.000000+5 2.513200+1 0.000000+0 1.030400-2 1.500000-1 0.000000+04531 2151  161\n"
    " 5.000000+0 0.000000+0          5          0        174         284531 2151  162\n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04531 2151  163\n"
    " 7.500000+0 2.170400+1 0.000000+0 8.898700-3 1.500000-1 0.000000+04531 2151  164\n"
    " 1.000000+1 2.079500+1 0.000000+0 8.526000-3 1.500000-1 0.000000+04531 2151  165\n"
    " 1.000000+2 2.163100+1 0.000000+0 8.868600-3 1.500000-1 0.000000+04531 2151  166\n"
    " 2.000000+2 2.176300+1 0.000000+0 8.922600-3 1.500000-1 0.000000+04531 2151  167\n"
    " 3.000000+2 2.182700+1 0.000000+0 8.949000-3 1.500000-1 0.000000+04531 2151  168\n"
    " 4.000000+2 2.187000+1 0.000000+0 8.966600-3 1.500000-1 0.000000+04531 2151  169\n"
    " 5.000000+2 2.189800+1 0.000000+0 8.978000-3 1.500000-1 0.000000+04531 2151  170\n"
    " 6.000000+2 2.192100+1 0.000000+0 8.987800-3 1.500000-1 0.000000+04531 2151  171\n"
    " 8.000000+2 2.195400+1 0.000000+0 9.001000-3 1.500000-1 0.000000+04531 2151  172\n"
    " 1.000000+3 2.197300+1 0.000000+0 9.009000-3 1.500000-1 0.000000+04531 2151  173\n"
    " 1.500000+3 2.201000+1 0.000000+0 9.024200-3 1.500000-1 0.000000+04531 2151  174\n"
    " 2.000000+3 2.203700+1 0.000000+0 9.035300-3 1.500000-1 0.000000+04531 2151  175\n"
    " 3.000000+3 2.207200+1 0.000000+0 9.049300-3 1.500000-1 0.000000+04531 2151  176\n"
    " 4.000000+3 2.199900+1 0.000000+0 9.019500-3 1.500000-1 0.000000+04531 2151  177\n"
    " 5.000000+3 2.202800+1 0.000000+0 9.031500-3 1.500000-1 0.000000+04531 2151  178\n"
    " 6.000000+3 2.205400+1 0.000000+0 9.042200-3 1.500000-1 0.000000+04531 2151  179\n"
    " 7.000000+3 2.206000+1 0.000000+0 9.044800-3 1.500000-1 0.000000+04531 2151  180\n"
    " 8.000000+3 2.207400+1 0.000000+0 9.050200-3 1.500000-1 0.000000+04531 2151  181\n"
    " 1.000000+4 2.209100+1 0.000000+0 9.057300-3 1.500000-1 0.000000+04531 2151  182\n"
    " 1.500000+4 2.208100+1 0.000000+0 9.053400-3 1.500000-1 0.000000+04531 2151  183\n"
    " 2.500000+4 2.204500+1 0.000000+0 9.038600-3 1.500000-1 0.000000+04531 2151  184\n"
    " 3.000000+4 2.196800+1 0.000000+0 9.006800-3 1.500000-1 0.000000+04531 2151  185\n"
    " 4.000000+4 2.178700+1 0.000000+0 8.932600-3 1.500000-1 0.000000+04531 2151  186\n"
    " 5.000000+4 2.158900+1 0.000000+0 8.851700-3 1.500000-1 0.000000+04531 2151  187\n"
    " 6.000000+4 2.137600+1 0.000000+0 8.764400-3 1.500000-1 0.000000+04531 2151  188\n"
    " 7.000000+4 2.116800+1 0.000000+0 8.679000-3 1.500000-1 0.000000+04531 2151  189\n"
    " 8.000000+4 2.096200+1 0.000000+0 8.594500-3 1.500000-1 0.000000+04531 2151  190\n"
    " 1.000000+5 2.056300+1 0.000000+0 8.430800-3 1.500000-1 0.000000+04531 2151  191\n"
    " 1.040050+2 0.000000+0          2          0          6          04531 2151  192\n"
    " 1.000000+0 0.000000+0          5          0        174         284531 2151  193\n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04531 2151  194\n"
    " 7.500000+0 7.958200+1 0.000000+0 4.456600-3 1.500000-1 0.000000+04531 2151  195\n"
    " 1.000000+1 7.624800+1 0.000000+0 4.269900-3 1.500000-1 0.000000+04531 2151  196\n"
    " 1.000000+2 7.931300+1 0.000000+0 4.441500-3 1.500000-1 0.000000+04531 2151  197\n"
    " 2.000000+2 7.979600+1 0.000000+0 4.468600-3 1.500000-1 0.000000+04531 2151  198\n"
    " 3.000000+2 8.003200+1 0.000000+0 4.481800-3 1.500000-1 0.000000+04531 2151  199\n"
    " 4.000000+2 8.019000+1 0.000000+0 4.490600-3 1.500000-1 0.000000+04531 2151  200\n"
    " 5.000000+2 8.029100+1 0.000000+0 4.496300-3 1.500000-1 0.000000+04531 2151  201\n"
    " 6.000000+2 8.037800+1 0.000000+0 4.501200-3 1.500000-1 0.000000+04531 2151  202\n"
    " 8.000000+2 8.049700+1 0.000000+0 4.507800-3 1.500000-1 0.000000+04531 2151  203\n"
    " 1.000000+3 8.056800+1 0.000000+0 4.511800-3 1.500000-1 0.000000+04531 2151  204\n"
    " 1.500000+3 8.070500+1 0.000000+0 4.519500-3 1.500000-1 0.000000+04531 2151  205\n"
    " 2.000000+3 8.080300+1 0.000000+0 4.525000-3 1.500000-1 0.000000+04531 2151  206\n"
    " 3.000000+3 8.092900+1 0.000000+0 4.532000-3 1.500000-1 0.000000+04531 2151  207\n"
    " 4.000000+3 8.066300+1 0.000000+0 4.517100-3 1.500000-1 0.000000+04531 2151  208\n"
    " 5.000000+3 8.076900+1 0.000000+0 4.523100-3 1.500000-1 0.000000+04531 2151  209\n"
    " 6.000000+3 8.086600+1 0.000000+0 4.528500-3 1.500000-1 0.000000+04531 2151  210\n"
    " 7.000000+3 8.088800+1 0.000000+0 4.529700-3 1.500000-1 0.000000+04531 2151  211\n"
    " 8.000000+3 8.093700+1 0.000000+0 4.532500-3 1.500000-1 0.000000+04531 2151  212\n"
    " 1.000000+4 8.100000+1 0.000000+0 4.536000-3 1.500000-1 0.000000+04531 2151  213\n"
    " 1.500000+4 8.096500+1 0.000000+0 4.534000-3 1.500000-1 0.000000+04531 2151  214\n"
    " 2.500000+4 8.083300+1 0.000000+0 4.526700-3 1.500000-1 0.000000+04531 2151  215\n"
    " 3.000000+4 8.054800+1 0.000000+0 4.510700-3 1.500000-1 0.000000+04531 2151  216\n"
    " 4.000000+4 7.988500+1 0.000000+0 4.473600-3 1.500000-1 0.000000+04531 2151  217\n"
    " 5.000000+4 7.916100+1 0.000000+0 4.433000-3 1.500000-1 0.000000+04531 2151  218\n"
    " 6.000000+4 7.838000+1 0.000000+0 4.389300-3 1.500000-1 0.000000+04531 2151  219\n"
    " 7.000000+4 7.761700+1 0.000000+0 4.346600-3 1.500000-1 0.000000+04531 2151  220\n"
    " 8.000000+4 7.686100+1 0.000000+0 4.304200-3 1.500000-1 0.000000+04531 2151  221\n"
    " 1.000000+5 7.539700+1 0.000000+0 4.222200-3 1.500000-1 0.000000+04531 2151  222\n"
    " 2.000000+0 0.000000+0          5          0        174         284531 2151  223\n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+04531 2151  224\n"
    " 7.500000+0 4.774900+1 0.000000+0 2.674000-3 1.500000-1 0.000000+04531 2151  225\n"
    " 1.000000+1 4.574900+1 0.000000+0 2.561900-3 1.500000-1 0.000000+04531 2151  226\n"
    " 1.000000+2 4.758800+1 0.000000+0 2.664900-3 1.500000-1 0.000000+04531 2151  227\n"
    " 2.000000+2 4.787800+1 0.000000+0 2.681100-3 1.500000-1 0.000000+04531 2151  228\n"
    " 3.000000+2 4.801900+1 0.000000+0 2.689100-3 1.500000-1 0.000000+04531 2151  229\n"
    " 4.000000+2 4.811400+1 0.000000+0 2.694400-3 1.500000-1 0.000000+04531 2151  230\n"
    " 5.000000+2 4.817500+1 0.000000+0 2.697800-3 1.500000-1 0.000000+04531 2151  231\n"
    " 6.000000+2 4.822700+1 0.000000+0 2.700700-3 1.500000-1 0.000000+04531 2151  232\n"
    " 8.000000+2 4.829800+1 0.000000+0 2.704700-3 1.500000-1 0.000000+04531 2151  233\n"
    " 1.000000+3 4.834100+1 0.000000+0 2.707100-3 1.500000-1 0.000000+04531 2151  234\n"
    " 1.500000+3 4.842300+1 0.000000+0 2.711700-3 1.500000-1 0.000000+04531 2151  235\n"
    " 2.000000+3 4.848200+1 0.000000+0 2.715000-3 1.500000-1 0.000000+04531 2151  236\n"
    " 3.000000+3 4.855700+1 0.000000+0 2.719200-3 1.500000-1 0.000000+04531 2151  237\n"
    " 4.000000+3 4.839800+1 0.000000+0 2.710300-3 1.500000-1 0.000000+04531 2151  238\n"
    " 5.000000+3 4.846200+1 0.000000+0 2.713900-3 1.500000-1 0.000000+04531 2151  239\n"
    " 6.000000+3 4.851900+1 0.000000+0 2.717100-3 1.500000-1 0.000000+04531 2151  240\n"
    " 7.000000+3 4.853300+1 0.000000+0 2.717800-3 1.500000-1 0.000000+04531 2151  241\n"
    " 8.000000+3 4.856200+1 0.000000+0 2.719500-3 1.500000-1 0.000000+04531 2151  242\n"
    " 1.000000+4 4.860000+1 0.000000+0 2.721600-3 1.500000-1 0.000000+04531 2151  243\n"
    " 1.500000+4 4.857900+1 0.000000+0 2.720400-3 1.500000-1 0.000000+04531 2151  244\n"
    " 2.500000+4 4.850000+1 0.000000+0 2.716000-3 1.500000-1 0.000000+04531 2151  245\n"
    " 3.000000+4 4.832900+1 0.000000+0 2.706400-3 1.500000-1 0.000000+04531 2151  246\n"
    " 4.000000+4 4.793100+1 0.000000+0 2.684100-3 1.500000-1 0.000000+04531 2151  247\n"
    " 5.000000+4 4.749700+1 0.000000+0 2.659800-3 1.500000-1 0.000000+04531 2151  248\n"
    " 6.000000+4 4.702800+1 0.000000+0 2.633600-3 1.500000-1 0.000000+04531 2151  249\n"
    " 7.000000+4 4.657000+1 0.000000+0 2.607900-3 1.500000-1 0.000000+04531 2151  250\n"
    " 8.000000+4 4.611700+1 0.000000+0 2.582500-3 1.500000-1 0.000000+04531 2151  251\n"
    " 1.000000+5 4.523800+1 0.000000+0 2.533300-3 1.500000-1 0.000000+04531 2151  252\n"
    " 3.000000+0 0.000000+0          5          0        174         284531 2151  253\n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+04531 2151  254\n"
    " 7.500000+0 3.410700+1 0.000000+0 1.910000-3 1.500000-1 0.000000+04531 2151  255\n"
    " 1.000000+1 3.267800+1 0.000000+0 1.830000-3 1.500000-1 0.000000+04531 2151  256\n"
    " 1.000000+2 3.399100+1 0.000000+0 1.903500-3 1.500000-1 0.000000+04531 2151  257\n"
    " 2.000000+2 3.419800+1 0.000000+0 1.915100-3 1.500000-1 0.000000+04531 2151  258\n"
    " 3.000000+2 3.429900+1 0.000000+0 1.920800-3 1.500000-1 0.000000+04531 2151  259\n"
    " 4.000000+2 3.436700+1 0.000000+0 1.924500-3 1.500000-1 0.000000+04531 2151  260\n"
    " 5.000000+2 3.441100+1 0.000000+0 1.927000-3 1.500000-1 0.000000+04531 2151  261\n"
    " 6.000000+2 3.444800+1 0.000000+0 1.929100-3 1.500000-1 0.000000+04531 2151  262\n"
    " 8.000000+2 3.449900+1 0.000000+0 1.931900-3 1.500000-1 0.000000+04531 2151  263\n"
    " 1.000000+3 3.452900+1 0.000000+0 1.933600-3 1.500000-1 0.000000+04531 2151  264\n"
    " 1.500000+3 3.458800+1 0.000000+0 1.936900-3 1.500000-1 0.000000+04531 2151  265\n"
    " 2.000000+3 3.463000+1 0.000000+0 1.939300-3 1.500000-1 0.000000+04531 2151  266\n"
    " 3.000000+3 3.468400+1 0.000000+0 1.942300-3 1.500000-1 0.000000+04531 2151  267\n"
    " 4.000000+3 3.457000+1 0.000000+0 1.935900-3 1.500000-1 0.000000+04531 2151  268\n"
    " 5.000000+3 3.461500+1 0.000000+0 1.938500-3 1.500000-1 0.000000+04531 2151  269\n"
    " 6.000000+3 3.465700+1 0.000000+0 1.940800-3 1.500000-1 0.000000+04531 2151  270\n"
    " 7.000000+3 3.466600+1 0.000000+0 1.941300-3 1.500000-1 0.000000+04531 2151  271\n"
    " 8.000000+3 3.468700+1 0.000000+0 1.942500-3 1.500000-1 0.000000+04531 2151  272\n"
    " 1.000000+4 3.471400+1 0.000000+0 1.944000-3 1.500000-1 0.000000+04531 2151  273\n"
    " 1.500000+4 3.469900+1 0.000000+0 1.943200-3 1.500000-1 0.000000+04531 2151  274\n"
    " 2.500000+4 3.464300+1 0.000000+0 1.940000-3 1.500000-1 0.000000+04531 2151  275\n"
    " 3.000000+4 3.452100+1 0.000000+0 1.933200-3 1.500000-1 0.000000+04531 2151  276\n"
    " 4.000000+4 3.423700+1 0.000000+0 1.917200-3 1.500000-1 0.000000+04531 2151  277\n"
    " 5.000000+4 3.392600+1 0.000000+0 1.899900-3 1.500000-1 0.000000+04531 2151  278\n"
    " 6.000000+4 3.359200+1 0.000000+0 1.881100-3 1.500000-1 0.000000+04531 2151  279\n"
    " 7.000000+4 3.326400+1 0.000000+0 1.862800-3 1.500000-1 0.000000+04531 2151  280\n"
    " 8.000000+4 3.294100+1 0.000000+0 1.844700-3 1.500000-1 0.000000+04531 2151  281\n"
    " 1.000000+5 3.231300+1 0.000000+0 1.809500-3 1.500000-1 0.000000+04531 2151  282\n"
    " 4.000000+0 0.000000+0          5          0        174         284531 2151  283\n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+04531 2151  284\n"
    " 7.500000+0 2.652700+1 0.000000+0 1.485500-3 1.500000-1 0.000000+04531 2151  285\n"
    " 1.000000+1 2.541600+1 0.000000+0 1.423300-3 1.500000-1 0.000000+04531 2151  286\n"
    " 1.000000+2 2.643800+1 0.000000+0 1.480500-3 1.500000-1 0.000000+04531 2151  287\n"
    " 2.000000+2 2.659900+1 0.000000+0 1.489500-3 1.500000-1 0.000000+04531 2151  288\n"
    " 3.000000+2 2.667700+1 0.000000+0 1.493900-3 1.500000-1 0.000000+04531 2151  289\n"
    " 4.000000+2 2.673000+1 0.000000+0 1.496900-3 1.500000-1 0.000000+04531 2151  290\n"
    " 5.000000+2 2.676400+1 0.000000+0 1.498800-3 1.500000-1 0.000000+04531 2151  291\n"
    " 6.000000+2 2.679300+1 0.000000+0 1.500400-3 1.500000-1 0.000000+04531 2151  292\n"
    " 8.000000+2 2.683200+1 0.000000+0 1.502600-3 1.500000-1 0.000000+04531 2151  293\n"
    " 1.000000+3 2.685600+1 0.000000+0 1.503900-3 1.500000-1 0.000000+04531 2151  294\n"
    " 1.500000+3 2.690200+1 0.000000+0 1.506500-3 1.500000-1 0.000000+04531 2151  295\n"
    " 2.000000+3 2.693400+1 0.000000+0 1.508300-3 1.500000-1 0.000000+04531 2151  296\n"
    " 3.000000+3 2.697600+1 0.000000+0 1.510700-3 1.500000-1 0.000000+04531 2151  297\n"
    " 4.000000+3 2.688800+1 0.000000+0 1.505700-3 1.500000-1 0.000000+04531 2151  298\n"
    " 5.000000+3 2.692300+1 0.000000+0 1.507700-3 1.500000-1 0.000000+04531 2151  299\n"
    " 6.000000+3 2.695500+1 0.000000+0 1.509500-3 1.500000-1 0.000000+04531 2151  300\n"
    " 7.000000+3 2.696300+1 0.000000+0 1.509900-3 1.500000-1 0.000000+04531 2151  301\n"
    " 8.000000+3 2.697900+1 0.000000+0 1.510800-3 1.500000-1 0.000000+04531 2151  302\n"
    " 1.000000+4 2.700000+1 0.000000+0 1.512000-3 1.500000-1 0.000000+04531 2151  303\n"
    " 1.500000+4 2.698800+1 0.000000+0 1.511300-3 1.500000-1 0.000000+04531 2151  304\n"
    " 2.500000+4 2.694400+1 0.000000+0 1.508900-3 1.500000-1 0.000000+04531 2151  305\n"
    " 3.000000+4 2.684900+1 0.000000+0 1.503600-3 1.500000-1 0.000000+04531 2151  306\n"
    " 4.000000+4 2.662800+1 0.000000+0 1.491200-3 1.500000-1 0.000000+04531 2151  307\n"
    " 5.000000+4 2.638700+1 0.000000+0 1.477700-3 1.500000-1 0.000000+04531 2151  308\n"
    " 6.000000+4 2.612700+1 0.000000+0 1.463100-3 1.500000-1 0.000000+04531 2151  309\n"
    " 7.000000+4 2.587200+1 0.000000+0 1.448900-3 1.500000-1 0.000000+04531 2151  310\n"
    " 8.000000+4 2.562000+1 0.000000+0 1.434700-3 1.500000-1 0.000000+04531 2151  311\n"
    " 1.000000+5 2.513200+1 0.000000+0 1.407400-3 1.500000-1 0.000000+04531 2151  312\n"
    " 5.000000+0 0.000000+0          5          0        174         284531 2151  313\n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+04531 2151  314\n"
    " 7.500000+0 2.170400+1 0.000000+0 1.215400-3 1.500000-1 0.000000+04531 2151  315\n"
    " 1.000000+1 2.079500+1 0.000000+0 1.164500-3 1.500000-1 0.000000+04531 2151  316\n"
    " 1.000000+2 2.163100+1 0.000000+0 1.211300-3 1.500000-1 0.000000+04531 2151  317\n"
    " 2.000000+2 2.176300+1 0.000000+0 1.218700-3 1.500000-1 0.000000+04531 2151  318\n"
    " 3.000000+2 2.182700+1 0.000000+0 1.222300-3 1.500000-1 0.000000+04531 2151  319\n"
    " 4.000000+2 2.187000+1 0.000000+0 1.224700-3 1.500000-1 0.000000+04531 2151  320\n"
    " 5.000000+2 2.189800+1 0.000000+0 1.226300-3 1.500000-1 0.000000+04531 2151  321\n"
    " 6.000000+2 2.192100+1 0.000000+0 1.227600-3 1.500000-1 0.000000+04531 2151  322\n"
    " 8.000000+2 2.195400+1 0.000000+0 1.229400-3 1.500000-1 0.000000+04531 2151  323\n"
    " 1.000000+3 2.197300+1 0.000000+0 1.230500-3 1.500000-1 0.000000+04531 2151  324\n"
    " 1.500000+3 2.201000+1 0.000000+0 1.232600-3 1.500000-1 0.000000+04531 2151  325\n"
    " 2.000000+3 2.203700+1 0.000000+0 1.234100-3 1.500000-1 0.000000+04531 2151  326\n"
    " 3.000000+3 2.207200+1 0.000000+0 1.236000-3 1.500000-1 0.000000+04531 2151  327\n"
    " 4.000000+3 2.199900+1 0.000000+0 1.231900-3 1.500000-1 0.000000+04531 2151  328\n"
    " 5.000000+3 2.202800+1 0.000000+0 1.233600-3 1.500000-1 0.000000+04531 2151  329\n"
    " 6.000000+3 2.205400+1 0.000000+0 1.235000-3 1.500000-1 0.000000+04531 2151  330\n"
    " 7.000000+3 2.206000+1 0.000000+0 1.235400-3 1.500000-1 0.000000+04531 2151  331\n"
    " 8.000000+3 2.207400+1 0.000000+0 1.236100-3 1.500000-1 0.000000+04531 2151  332\n"
    " 1.000000+4 2.209100+1 0.000000+0 1.237100-3 1.500000-1 0.000000+04531 2151  333\n"
    " 1.500000+4 2.208100+1 0.000000+0 1.236600-3 1.500000-1 0.000000+04531 2151  334\n"
    " 2.500000+4 2.204500+1 0.000000+0 1.234500-3 1.500000-1 0.000000+04531 2151  335\n"
    " 3.000000+4 2.196800+1 0.000000+0 1.230200-3 1.500000-1 0.000000+04531 2151  336\n"
    " 4.000000+4 2.178700+1 0.000000+0 1.220100-3 1.500000-1 0.000000+04531 2151  337\n"
    " 5.000000+4 2.158900+1 0.000000+0 1.209000-3 1.500000-1 0.000000+04531 2151  338\n"
    " 6.000000+4 2.137600+1 0.000000+0 1.197100-3 1.500000-1 0.000000+04531 2151  339\n"
    " 7.000000+4 2.116800+1 0.000000+0 1.185400-3 1.500000-1 0.000000+04531 2151  340\n"
    " 8.000000+4 2.096200+1 0.000000+0 1.173900-3 1.500000-1 0.000000+04531 2151  341\n"
    " 1.000000+5 2.056300+1 0.000000+0 1.151500-3 1.500000-1 0.000000+04531 2151  342\n"
    " 6.000000+0 0.000000+0          5          0        174         284531 2151  343\n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04531 2151  344\n"
    " 7.500000+0 1.836500+1 0.000000+0 1.028400-3 1.500000-1 0.000000+04531 2151  345\n"
    " 1.000000+1 1.759600+1 0.000000+0 9.853600-4 1.500000-1 0.000000+04531 2151  346\n"
    " 1.000000+2 1.830300+1 0.000000+0 1.025000-3 1.500000-1 0.000000+04531 2151  347\n"
    " 2.000000+2 1.841400+1 0.000000+0 1.031200-3 1.500000-1 0.000000+04531 2151  348\n"
    " 3.000000+2 1.846900+1 0.000000+0 1.034300-3 1.500000-1 0.000000+04531 2151  349\n"
    " 4.000000+2 1.850500+1 0.000000+0 1.036300-3 1.500000-1 0.000000+04531 2151  350\n"
    " 5.000000+2 1.852900+1 0.000000+0 1.037600-3 1.500000-1 0.000000+04531 2151  351\n"
    " 6.000000+2 1.854900+1 0.000000+0 1.038700-3 1.500000-1 0.000000+04531 2151  352\n"
    " 8.000000+2 1.857600+1 0.000000+0 1.040300-3 1.500000-1 0.000000+04531 2151  353\n"
    " 1.000000+3 1.859300+1 0.000000+0 1.041200-3 1.500000-1 0.000000+04531 2151  354\n"
    " 1.500000+3 1.862400+1 0.000000+0 1.043000-3 1.500000-1 0.000000+04531 2151  355\n"
    " 2.000000+3 1.864700+1 0.000000+0 1.044200-3 1.500000-1 0.000000+04531 2151  356\n"
    " 3.000000+3 1.867600+1 0.000000+0 1.045900-3 1.500000-1 0.000000+04531 2151  357\n"
    " 4.000000+3 1.861400+1 0.000000+0 1.042400-3 1.500000-1 0.000000+04531 2151  358\n"
    " 5.000000+3 1.863900+1 0.000000+0 1.043800-3 1.500000-1 0.000000+04531 2151  359\n"
    " 6.000000+3 1.866100+1 0.000000+0 1.045000-3 1.500000-1 0.000000+04531 2151  360\n"
    " 7.000000+3 1.866700+1 0.000000+0 1.045300-3 1.500000-1 0.000000+04531 2151  361\n"
    " 8.000000+3 1.867800+1 0.000000+0 1.046000-3 1.500000-1 0.000000+04531 2151  362\n"
    " 1.000000+4 1.869200+1 0.000000+0 1.046800-3 1.500000-1 0.000000+04531 2151  363\n"
    " 1.500000+4 1.868400+1 0.000000+0 1.046300-3 1.500000-1 0.000000+04531 2151  364\n"
    " 2.500000+4 1.865400+1 0.000000+0 1.044600-3 1.500000-1 0.000000+04531 2151  365\n"
    " 3.000000+4 1.858800+1 0.000000+0 1.040900-3 1.500000-1 0.000000+04531 2151  366\n"
    " 4.000000+4 1.843500+1 0.000000+0 1.032400-3 1.500000-1 0.000000+04531 2151  367\n"
    " 5.000000+4 1.826800+1 0.000000+0 1.023000-3 1.500000-1 0.000000+04531 2151  368\n"
    " 6.000000+4 1.808800+1 0.000000+0 1.012900-3 1.500000-1 0.000000+04531 2151  369\n"
    " 7.000000+4 1.791200+1 0.000000+0 1.003100-3 1.500000-1 0.000000+04531 2151  370\n"
    " 8.000000+4 1.773700+1 0.000000+0 9.932800-4 1.500000-1 0.000000+04531 2151  371\n"
    " 1.000000+5 1.739900+1 0.000000+0 9.743600-4 1.500000-1 0.000000+04531 2151  372\n";
}
