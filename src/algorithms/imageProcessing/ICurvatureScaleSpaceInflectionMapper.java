package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.List;

/**
 *
 * @author nichole
 */
public interface ICurvatureScaleSpaceInflectionMapper {

    TransformationParameters createEuclideanTransformation();

    PairInt[] getMatchedEdgesIndexes();

    PairIntArray getMatchedXY1();

    float[] getMatchedXY1Weights();

    PairIntArray getMatchedXY2();

    float[] getMatchedXY2Weights();

    void initialize();

    void setToRefineTransformations();

    void useDebugMode();

    void useLineDrawingLineMode();

    void useOutdoorMode();
    
    public List<List<CurvatureScaleSpaceContour>> getContours1();
    
    public List<List<CurvatureScaleSpaceContour>> getContours2();
    
}
