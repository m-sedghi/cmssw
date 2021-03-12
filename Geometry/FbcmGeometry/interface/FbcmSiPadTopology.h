///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#ifndef Geometry_FbcmGeometry_FbcmSiPadTopology_H
#define Geometry_FbcmGeometry_FbcmSiPadTopology_H

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelChannelIdentifier.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class FbcmSiPadTopology final : public PixelTopology {
public:
  // Constructor, initilize
  FbcmSiPadTopology(float pitchx,
                    float pitchy )
      : m_pitchx(pitchx),
        m_pitchy(pitchy),
        m_nrows(1),
        m_ncols(1),
        m_ROWS_PER_ROC(1),  // Num of Rows per ROC 
        m_COLS_PER_ROC(1),  // Num of Cols per ROC 
        m_ROCS_X(1),        //  Number of ROCS X 
        m_ROCS_Y(1)        //  Number of ROCS Y 
		{
    // Calculate the edge of the active sensor with respect to the center,
    // that is simply the half-size.
    // Take into account large pixels
    //m_xoffset = -(m_nrows + BIG_PIX_PER_ROC_X * m_nrows / ROWS_PER_ROC) / 2. * m_pitchx;
    //m_yoffset = -(m_ncols + BIG_PIX_PER_ROC_Y * m_ncols / COLS_PER_ROC) / 2. * m_pitchy;
	
	m_xoffset = -m_nrows / 2. * m_pitchx;
    m_yoffset = -m_ncols / 2. * m_pitchy;

//std::cout << "m_xoffset:"<< m_xoffset << ", m_yoffset:"<< m_yoffset<<"\n"; 

    LogDebug("FbcmSiPadTopology") 		 << "nrows " << m_nrows 
										 << ", ncols " << m_ncols 
										 << ", pitchx " << m_pitchx
                                         << ", pitchy " << m_pitchy 
										 << ", xoffset " << m_xoffset 
										 << ", yoffset " << m_yoffset 
										 << ", ROWS_PER_ROC " << m_ROWS_PER_ROC 
										 << ", COLS_PER_ROC " << m_COLS_PER_ROC 
										 << ", ROCS_X " << m_ROCS_X
                                         << ", ROCS_Y " << m_ROCS_Y 
										 << "\nNROWS " << m_ROWS_PER_ROC * m_ROCS_X
                                         << ", NCOL " << m_COLS_PER_ROC * m_ROCS_Y;
  }

  // Topology interface, go from Masurement to Local corrdinates
  // SiPad coordinates (mp) -> cm (LocalPoint)
  LocalPoint localPosition(const MeasurementPoint& mp) const override;

  // Transform LocalPoint to Measurement. Call SiPad().
  MeasurementPoint measurementPosition(const LocalPoint& lp) const override {
    std::pair<float, float> p = SiPad(lp);
    return MeasurementPoint(p.first, p.second);
  }

  // PixelTopology interface.
  // Transform LocalPoint in cm to measurement in pitch units.
  std::pair<float, float> pixel(const LocalPoint& p) const override;
  std::pair<float, float> SiPad(const LocalPoint& p) const;

  // Errors
  // Error in local (cm) from the masurement errors
  LocalError localError(const MeasurementPoint&, const MeasurementError&) const override;
  // Errors in pitch units from localpoint error (in cm)
  MeasurementError measurementError(const LocalPoint&, const LocalError&) const override;

  //-------------------------------------------------------------
  // Transform LocalPoint to channel. Call SiPad()
  //
  int channel(const LocalPoint& lp) const override {
    std::pair<float, float> p = SiPad(lp);
    return PixelChannelIdentifier::pixelToChannel(int(p.first), int(p.second));
  }

  //-------------------------------------------------------------
  // Transform measurement to local coordinates individually in each dimension
  //
  float localX(const float mpX) const override;
  float localY(const float mpY) const override;

  //-------------------------------------------------------------
  // There is No BIG SiPad !!
  bool isItBigPixelInX(const int ixbin) const override { return false; }
  bool isItBigPixelInY(const int iybin) const override { return false; }
  bool containsBigPixelInX(int ixmin, int ixmax) const override { return false; }
  bool containsBigPixelInY(int iymin, int iymax) const override { return false; }
  

  //-------------------------------------------------------------
  // Check whether the SiPad is at the edge of the SilionDie
  //
  bool isItEdgePixelInX(int ixbin) const override { return ((ixbin == 0) | (ixbin == (m_nrows - 1))); }
  bool isItEdgePixelInY(int iybin) const override { return ((iybin == 0) | (iybin == (m_ncols - 1))); }
  bool isItEdgePixel(int ixbin, int iybin) const override { return (isItEdgePixelInX(ixbin) | isItEdgePixelInY(iybin));  }

  //------------------------------------------------------------------
  // Return pitch
  std::pair<float, float> pitch() const override { return std::pair<float, float>(float(m_pitchx), float(m_pitchy)); }
  // Return number of rows
  int nrows() const override { return (m_nrows); }
  // Return number of cols
  int ncolumns() const override { return (m_ncols); }
  // mlw Return number of ROCS Y
  int rocsY() const override { return m_ROCS_Y; }
  // mlw Return number of ROCS X
  int rocsX() const override { return m_ROCS_X; }
  // mlw Return number of rows per roc
  int rowsperroc() const override { return m_ROWS_PER_ROC; }
  // mlw Return number of cols per roc
  int colsperroc() const override { return m_COLS_PER_ROC; }
  float xoffset() const { return m_xoffset; }
  float yoffset() const { return m_yoffset; }

private:
  float m_pitchx;
  float m_pitchy;
  float m_xoffset;
  float m_yoffset;
  int m_nrows;
  int m_ncols;
  int m_ROWS_PER_ROC;
  int m_COLS_PER_ROC;
  int m_ROCS_X;
  int m_ROCS_Y;
};

#endif
