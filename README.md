_**This is the supporting information package of the codes and data for the paper "Multi-stable Spatial linkages".**_

**Paper title:**   
Multi-stable Spatial linkages

**Authors:**   
Tong Zhou[1], Zhuangzhi Miao[1], Chong Huang[1], Yang Li*[1,2] (yang.li@whu.edu.cn)

**Affiliations:**   
[1] The Institute of Technological Sciences, Wuhan University
[2] Wuhan University Shenzhen Research Institute

**Abstract:**
Multi-stable structures can be reconfigured with fewer, lightweight, and less accurate actuators. This is because the attraction domain in the multi-stable energy landscape provides both reconfiguration guidance and shape accuracy. Additionally, such structures can generate impulsive motion due to structural instability to mimic Venus flytraps and jumping insects. Most multi-stable structures are designed based on the geometric symmetry (bi-stable unit cells), or direct assembling of bi-stable unit cells (multi-stable assembly). Most multi-stable unit cells are planar structures, while spatial linkages can generate complex 3D motion and hold a more promising potential for application. This paper proposes a generalized approach to design a type of multi-stable spatial linkages (MSSLs) with multiple prescriptible configurations, which are structurally compatible, and naturally stable at these states. It reveals that all over-constrained mechanisms can be transformed into multi-stable structures with the same design method. Single-loop bi-stable 4R and tri-stable 6R spatial linkages modules with non-symmetric stable states, which are respectively transformed from Bennett and Bricard linkages, are designed to illustrate the intrinsic superiority over the ordinary methods, and then assembled deliberately into multi-loop linkages for practical aims. Three bi-stable multi-loop assemblies—a reconfigurable tube, a bio-inspired impulsive gripper, and a swimming robot—are experimentally demonstrated. This design method paves the way for morphing platforms with lightweight actuation, high shape accuracy, high stiffness, and prescribed impulsive 3D motion.

**Contents in this package:**
1. **Data S1**. Codes package to generate the single-loop cases of 4R and 6R MSSLs  
_The codes to generate the single-loop cases are provided. Run file named “MSSL_Bistable4R.m” and “MSSL_Tristable6R.m” to calculate respondingly._

2. **Data S2**. Stiffness analysis with MATLAB language  
_This code presents a detailed analysis for the stiffness in a 4R MSSL. Adopting this algorithm in their own design is highly recommended._

3. **Data S3**. Derivation of configuration constraints with Wolfram language
_This code presents a detailed derivation for the closure conditions in a 4R MSSL. There is a relatively clean and explicit formation for the equations and the incompatibility._

4. **Data S4**. G-J-K algorithm with MATLAB language  
_It is acknowledged that the code is created by Matthew Sheen (2023). (Fast 3D Collision Detection -- GJK algorithm (https://github.com/mws262/MATLAB-GJK-Collision-Detection ), GitHub. Retrieved November 2, 2023.)_

5. **Data S5**. Unwrapping algorithm of 4R and 6R MSSLs  
_This code gives an introductory example of unwrapping a 4R and 6R MSSL to a flat pattern for fast fabrication._

6. **Data S6**. Pattern of 4R and 6R MSSLs with different design targets  
_Apart from the classic 4R MSSL design in Fig. S7, there are other different design of 4R and 6R MSSLs provided in the dataset._

7. **Data S7**. Detailed technical data sheet of carbon fibers and TPU  
_Two data sheets of carbon fiber (PolyMide™ PA12-CF) and TPU (PolyFlex™ TPU95) are provided by technical support of UltiMaker._

8. **Data S8**. Data package of design results.  
_It contains the data in the parametric study (“GeneratedResult_ID1000.csv”), the simulation motion data of 4R MSSL (“MATX_Points_Motion_Bistable4R.csv”. Every 6 rows is a group of snapshots of motion), and the designs of 4R and 6R MSSL. (Files is with the name starting by “MSSL_”. Variable “MATX_ParaStruc” stores the D-H parameters. Variable “MATX_PointsS” stores the points at the stable states. Variable “MATX_SlipDist” stores the sliding distances.)_

9. **Data S9**. Codes package of figure in the paper  
_It contains the codes for the figure in the main text to review._
