#include <iostream>
#include <fstream>
#include <sstream> 

#include "edgeconfig.hpp"
#include "util.hpp"

#include <QtXml>
#include <QFile>
#include <QString>

#include <tiffio.h>

using namespace std;

template <typename T> T fromString(const QString& s, T &defval) {  return util::fromString(s.toStdString(), defval); }

class EdgeConfigHandler : public QXmlDefaultHandler {
  edgeconfig &conf;
public:
  EdgeConfigHandler(edgeconfig &conf_) : conf(conf_) { }
  bool startElement ( const QString & /*namespaceURI*/, const QString & /*localName*/, 
		      const QString & qName, const QXmlAttributes & attr ) {
    if(qName == "system") {
      conf.threads = fromString(attr.value("threads"), conf.threads);
      conf.hs_threads = fromString(attr.value("hs_threads"), conf.hs_threads);
      conf.analysis_id = fromString(attr.value("analysis_id"), conf.analysis_id);
      conf.run_refinement = attr.value("run_refinement").toInt() != 0;
      conf.run_mesh = attr.value("run_mesh").toInt() != 0;
      conf.keep_edge_channel_only = attr.value("keep_edge_channel_only").toInt() != 0;
      conf.keep_processing = attr.value("keep_processing").toInt() != 0;
      conf.edge_channel = fromString(attr.value("edge_channel"), conf.edge_channel);
      conf.nuc_channel = fromString(attr.value("nuc_channel"), conf.nuc_channel);
      conf.use_GL_buffers = attr.value("use_GL_buffers").toInt() != 0;
    }
    else if(qName == "tiff") {
      conf.frames = fromString(attr.value("frames"), conf.frames);
      conf.slices = fromString(attr.value("slices"), conf.slices);
      conf.channels = fromString(attr.value("channels"), conf.channels);
      conf.override_tiff = fromString(attr.value("override_tiff"), conf.override_tiff);
      conf.startframe = fromString(attr.value("startframe"), conf.startframe);
      conf.endframe = fromString(attr.value("endframe"), conf.endframe);
      conf.stepframe = fromString(attr.value("stepframe"), conf.stepframe);
      conf.voxelXY = fromString(attr.value("voxelXY"), conf.voxelXY);
      conf.voxelX = fromString(attr.value("voxelX"), conf.voxelXY); // Use x,y voxel as default (for old files).
      conf.voxelY = fromString(attr.value("voxelY"), conf.voxelXY);
      conf.voxelXY = min(conf.voxelX, conf.voxelY);
      cout << "voxelX=" << conf.voxelX << endl;
      cout << "voxelX=" << attr.value("voxelX").toStdString() << endl;
      cout << "voxelY=" << conf.voxelY << endl;
      conf.voxelZ = fromString(attr.value("voxelZ"), conf.voxelZ);
      conf.scalef = fromString(attr.value("scalef"), conf.scalef);
      conf.scaleToZ = fromString(attr.value("scaleToZ"), conf.scaleToZ);
      conf.time_step = fromString(attr.value("time_step"), conf.time_step);
    }
    else if(qName == "enhance") {
      conf.hist_nx = fromString(attr.value("hist_nx"), conf.hist_nx);
      conf.hist_ny = fromString(attr.value("hist_ny"), conf.hist_ny);
      conf.hist_nz = fromString(attr.value("hist_nz"), conf.hist_nz);
      conf.fmin = fromString(attr.value("fmin"), conf.fmin);
      conf.fmax = fromString(attr.value("fmax"), conf.fmax);
      // Membrane enhancement
      conf.r1 = fromString(attr.value("r1"), conf.r1);
      conf.r2 = fromString(attr.value("r2"), conf.r2);
      conf.alpha1 = fromString(attr.value("alpha1"), conf.alpha1);
      conf.alpha2 = fromString(attr.value("alpha2"), conf.alpha2);
      // Nucleus enhancement
      conf.r1_nuc = fromString(attr.value("r1_nuc"), conf.r1_nuc);
      conf.r2_nuc = fromString(attr.value("r2_nuc"), conf.r2_nuc);
      conf.alpha1_nuc = fromString(attr.value("alpha1_nuc"), conf.alpha1_nuc);
      conf.alpha2_nuc = fromString(attr.value("alpha2_nuc"), conf.alpha2_nuc);
      conf.r_nonedge = fromString(attr.value("r_nonedge"), conf.r_nonedge);
      conf.repair_nucs = attr.value("repair_nucs").toInt() != 0;
    }
    else if(qName == "dog") {
      conf.sigma_low = fromString(attr.value("sigma_low"), conf.sigma_low);
      conf.sigma_nuc = fromString(attr.value("sigma_nuc"), conf.sigma_nuc);

      conf.thresh_low = fromString(attr.value("thresh_low"), conf.thresh_low);
      conf.thresh_nuc = fromString(attr.value("thresh_nuc"), conf.thresh_nuc);

    }
    else if(qName == "mesh") {
      conf.mcdim        = fromString(attr.value("mcdim"), conf.mcdim);

      conf.min_nuc_vol = fromString(attr.value("min_nuc_vol"), conf.min_nuc_vol);
      conf.max_nuc_vol = fromString(attr.value("max_nuc_vol"), conf.max_nuc_vol);
      conf.nuc_split_dist = fromString(attr.value("nuc_split_dist"), conf.nuc_split_dist);

      // Segmentation params.
      conf.max_hole_rad = fromString(attr.value("max_hole_rad"), conf.max_hole_rad);
      conf.min_comp_vol = fromString(attr.value("min_comp_vol"), conf.min_comp_vol);
      conf.max_comp_vol = fromString(attr.value("max_comp_vol"), conf.max_comp_vol);
      conf.internal_limit = fromString(attr.value("internal_limit"), conf.internal_limit);
      conf.noise_comp_vol = fromString(attr.value("noise_comp_vol"), conf.noise_comp_vol);

      conf.refine_sigma = fromString(attr.value("refine_sigma"), conf.refine_sigma);
      conf.refine_dist  = fromString(attr.value("refine_dist"), conf.refine_dist);
      conf.refine_alpha = fromString(attr.value("refine_alpha"), conf.refine_alpha);
      conf.refine_iter  = fromString(attr.value("refine_iter"), conf.refine_iter);
      cout << "conf.refine_iter=" << conf.refine_iter << endl;
      conf.refine_stepsize = fromString(attr.value("refine_stepsize"), conf.refine_stepsize);
      conf.neighbor_alpha = fromString(attr.value("neighbor_alpha"), conf.neighbor_alpha);
    }
    else if(qName == "NM2010") {
      conf.run_NM2010 = attr.value("run_NM2010").toInt() != 0;
      conf.use_ASF = attr.value("use_ASF").toInt() != 0;
      conf.maxASF_rad = fromString(attr.value("maxASF_rad"), conf.maxASF_rad);	
      conf.filter_sigma = fromString(attr.value("filter_sigma"), conf.filter_sigma);	
      conf.bgseed_thresh = fromString(attr.value("bgseed_thresh"), conf.bgseed_thresh);
      conf.h_seed = fromString(attr.value("h_seed"), conf.h_seed);
      conf.error_vox_thresh = fromString(attr.value("error_vox_thresh"), conf.error_vox_thresh);
    }
    else if(qName == "boundary") {
      boundaryconf_t b_conf;
      b_conf.sigma_high = fromString(attr.value("sigma_high"), b_conf.sigma_high);
      b_conf.thresh_high = fromString(attr.value("thresh_high"), b_conf.thresh_high);
      b_conf.boundary_fmin = fromString(attr.value("boundary_fmin"), b_conf.boundary_fmin);
      b_conf.boundary_mincomp = fromString(attr.value("boundary_mincomp"), b_conf.boundary_mincomp);
      b_conf.r_repair = fromString(attr.value("r_repair"), b_conf.r_repair);
      b_conf.r2_repair = fromString(attr.value("r2_repair"), b_conf.r2_repair);
      b_conf.startframe = fromString(attr.value("startframe"), b_conf.startframe);
      b_conf.endframe = fromString(attr.value("endframe"), b_conf.endframe);
      if(b_conf.startframe == 0 && b_conf.endframe == 0) conf.boundaryconf.at(0) = b_conf;
      else conf.boundaryconf.push_back(b_conf);
    }
    else if(qName == "ACME") {
      conf.ACME_average_radius = fromString(attr.value("average_radius"), conf.ACME_average_radius);
      conf.ACME_plate_measure = fromString(attr.value("plate_measure"), conf.ACME_plate_measure);
      conf.ACME_voting_field = fromString(attr.value("voting_field"), conf.ACME_voting_field);
      conf.ACME_membrane_segmentation = fromString(attr.value("membrane_segmentation"), conf.ACME_membrane_segmentation);
    }
    return true;
  }
  bool endElement ( const QString & /*namespaceURI*/, const QString & /*localName*/, const QString & /*qName*/ ) { return true; }
  bool characters ( const QString & /*ch*/ ) { return true; }
};


bool edgeconfig::load_xml(string filename) {
  QXmlSimpleReader xmlReader;
  QFile *file = new QFile(filename.c_str());
  QXmlInputSource *source = new QXmlInputSource(file);
  EdgeConfigHandler *handler = new EdgeConfigHandler(*this);
  xmlReader.setContentHandler(handler);
  bool ret = xmlReader.parse(source);
  delete source; delete handler; delete file; 
  return ret;
}

// Default handling of the output of a binary paramter. 
void B(ofstream &config, string key, bool value) { 
  if(value) config << key << "=\"1\" ";  else config << key << "=\"0\" "; 
}

template <typename T> 
void V(ofstream&config, string key, T value) {
  config << key << "=\"" << value << "\" ";
}

bool edgeconfig::save_xml(string filename) {
  ofstream conf(filename.c_str());
  if(!conf) return false;

  conf << "<edge>" << endl;

  conf << "<system ";
  conf << "threads=\"" << threads << "\" ";
  V(conf, "hs_threads", hs_threads);
  conf << "analysis_id=\"" << analysis_id << "\" ";
  B(conf, "run_refinement", run_refinement);
  B(conf, "run_mesh", run_mesh);
  B(conf, "keep_edge_channel_only", keep_edge_channel_only);
  conf << "edge_channel=\"" << edge_channel << "\" ";
  conf << "nuc_channel=\"" << nuc_channel << "\" ";
  B(conf, "keep_processing", keep_processing);
  B(conf, "use_GL_buffers", use_GL_buffers);
  conf << "/>" << endl << endl;
  
  // Write XML file. 
  conf << "<tiff ";
  conf << "frames=\"" << frames << "\" ";
  conf << "slices=\"" << slices << "\" ";
  conf << "channels=\"" << channels << "\" ";
  conf << "endframe=\"" << endframe << "\" ";
  B(conf, "override_tiff", override_tiff);
  conf << "startframe=\"" << startframe << "\" ";
  conf << "stepframe=\"" << stepframe << "\" ";
  conf << "voxelXY=\"" << voxelXY << "\" ";
  conf << "voxelX=\"" << voxelX << "\" ";
  conf << "voxelY=\"" << voxelY << "\" ";
  conf << "voxelZ=\"" << voxelZ << "\" ";
  conf << "scalef=\"" << scalef << "\" ";
  B(conf, "scaleToZ", scaleToZ);
  conf << "time_step=\"" << time_step << "\" ";
  conf << "/>" << endl << endl;;
  
  conf << "<enhance ";
  conf << "hist_nx=\"" << hist_nx << "\" ";
  conf << "hist_ny=\"" << hist_ny << "\" ";
  conf << "hist_nz=\"" << hist_nz << "\" ";
  conf << "fmin=\"" << fmin << "\" ";
  conf << "fmax=\"" << fmax << "\" ";
  conf << "r1=\"" << r1 << "\" ";
  conf << "r2=\"" << r2 << "\" ";
  conf << "alpha1=\"" << alpha1 << "\" ";
  conf << "alpha2=\"" << alpha2 << "\" ";
  conf << "r1_nuc=\"" << r1_nuc << "\" ";
  conf << "r2_nuc=\"" << r2_nuc << "\" ";
  conf << "alpha1_nuc=\"" << alpha1_nuc << "\" ";
  conf << "alpha2_nuc=\"" << alpha2_nuc << "\" ";
  conf << "r_nonedge=\"" << r_nonedge << "\" ";
  B(conf, "repair_nucs", repair_nucs);
  conf << "/>" << endl << endl;


  conf << "<dog ";
  conf << "sigma_low=\"" << sigma_low << "\" ";
  conf << "sigma_nuc=\"" << sigma_nuc << "\" ";
  conf << "thresh_low=\"" << thresh_low << "\" ";
  conf << "thresh_nuc=\"" << thresh_nuc << "\" ";
  conf << "/>" << endl << endl;

  for(size_t b = 0; b < boundaryconf.size(); b++) {
    boundaryconf_t &bconf = boundaryconf[b];
    conf << "<boundary ";
    V(conf, "startframe" ,bconf.startframe);
    V(conf, "endframe" ,bconf.endframe);
    V(conf, "boundary_mincomp", bconf.boundary_mincomp);
    V(conf, "boundary_fmin", bconf.boundary_fmin);
    V(conf, "thresh_high", bconf.thresh_high);
    V(conf, "sigma_high", bconf.sigma_high);
    V(conf, "r_repair", bconf.r_repair);
    V(conf, "r2_repair", bconf.r2_repair);
    conf << "/>" << endl << endl;
  }

  conf << "<mesh ";
  conf << "mcdim=\"" << mcdim << "\" ";
  conf << "min_nuc_vol=\"" << min_nuc_vol << "\" ";
  conf << "max_nuc_vol=\"" << max_nuc_vol << "\" ";
  conf << "nuc_split_dist=\"" << nuc_split_dist << "\" ";
  conf << "max_hole_rad=\"" << max_hole_rad << "\" ";
  conf << "min_comp_vol=\"" << min_comp_vol << "\" ";
  conf << "max_comp_vol=\"" << max_comp_vol << "\" ";
  conf << "internal_limit=\"" << internal_limit << "\" ";
  conf << "refine_sigma=\"" << refine_sigma << "\" ";
  conf << "refine_dist=\"" << refine_dist << "\" ";
  conf << "refine_alpha=\"" << refine_alpha << "\" ";
  conf << "refine_iter=\"" << refine_iter << "\" ";
  conf << "refine_stepsize=\"" << refine_stepsize << "\" ";
  conf << "noise_comp_vol=\"" << noise_comp_vol << "\" ";
  V(conf, "neighbor_alpha", neighbor_alpha);
  conf << "/>" << endl << endl;

  conf << "<NM2010 ";
  B(conf, "run_NM2010", run_NM2010);
  B(conf, "use_ASF", use_ASF);
  V(conf, "maxASF_rad", maxASF_rad);
  V(conf, "filter_sigma", filter_sigma);
  V(conf, "bgseed_thresh", bgseed_thresh);
  V(conf, "h_seed", h_seed);
  V(conf, "error_vox_thresh", error_vox_thresh);
  conf << "/>" << endl << endl;

  conf << "<ACME ";
  V(conf, "average_radius", ACME_average_radius);
  V(conf, "plate_measure", ACME_plate_measure);
  V(conf, "voting_field", ACME_voting_field);
  V(conf, "membrane_segmentation", ACME_membrane_segmentation);
  conf << "/>" << endl << endl;


  conf << "</edge>" << endl;

  return true;
}


bool edgeconfig::load_tiff(string tiffstack) {
  TIFFSetWarningHandler(NULL); TIFFSetErrorHandler(NULL);
  // Determine image stack depth.
  TIFF *tif = TIFFOpen(tiffstack.c_str(), "r");
  if(!tif) return false;
  // Set default values.
  slices = 0; frames = 1; channels = 1; 
  char *raw_description = NULL;
  TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &raw_description); 
  string descr = raw_description;
  bool has_slices = false;
  size_t i = 0; 
  while(i < descr.length()) {
    string key, value;
    while(i < descr.length() && descr[i] != '=') {
      if(descr[i] != ' ') key += descr[i]; i++;
    }
    i++;
    while(i < descr.length() && descr[i] != '\n') {
      if(descr[i] != ' ') value += descr[i]; i++;
    }
    if(key == "frames") { frames = util::fromString<int>(value); }
    if(key == "slices") { slices = util::fromString<int>(value); has_slices = true; }
    if(key == "channels") { channels = util::fromString<int>(value); }
    i++;
  }
  TIFFClose(tif);

  if(!has_slices) { // No slices, just one frame.
    tif = TIFFOpen(tiffstack.c_str(), "r");
    slices = 0;
    do { slices++; } while(TIFFReadDirectory(tif));
    TIFFClose(tif);
    slices /= channels;
  }
  cout << "slices=" << slices << " frames=" << frames << " channels=" << channels << endl;
  if(startframe < 0 || startframe > frames-1) startframe = 0;
  if(endframe < 0 || endframe > frames-1) endframe = frames-1;
  return true;
}

