--- src/load_Analysis.cc.txt	2015-02-02 12:32:50.000000000 +0100
+++ src/load_Analysis.cc	2015-02-02 12:33:03.000000000 +0100
@@ -244,7 +244,7 @@
     // Set up DijetAnalyzer
     AnalyzerDijet* analyzerDijet(0);
     if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::dijet) != v_analysisMode.end()){
-        analyzerDijet = new AnalyzerDijet(Mva2dWeightsFILE, "correct_step7_cate0_cate1_cate2_d144", "", {}, {"7"}, jetCategories, false, true);
+        analyzerDijet = new AnalyzerDijet(Mva2dWeightsFILE, "correct_step7_cate0_cate1_cate2_d144", "", {"0b"}, {"7"}, jetCategories, false, true);
         v_analyzer.push_back(analyzerDijet);
     }
     
