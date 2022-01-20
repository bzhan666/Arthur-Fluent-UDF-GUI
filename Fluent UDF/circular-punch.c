/********************************************************************
    UDF to demonstrate deforming mesh / grid motion
    Written by J Kidger

    A circular punch is pushed into a square domain

    Square is assumed to be 2m x 2m, centred at origin

   Preview with delta-T of 0.01s
   UDF演示变形的网格/网格运动
     由J Kidger撰写

     将圆形冲头推入方形区域

     假定正方形为2m x 2m，以原点为中心

    以0.01s的delta-T预览

********************************************************************/
#include "udf.h"

#define  circ_rad     0.75      /* radius of circular punch    圆冲头半径 */
#define  punch_y_vel  -1.0    /* velocity of punch moving down 冲头向下移动的速度*/


DEFINE_GRID_MOTION(jk_punch, domain, dt, time, dtime)
{
  Thread *tf = DT_THREAD (dt);
  face_t f;
  Node *node_p;
  real x0,x,y0,y1,y2;
  int n;

/* Compute start co-ordinates of punch, which starts just above square
计算冲头的开始坐标，该坐标刚好在正方形上方*/
   x0=0.0 ;
   y0=1.0 + circ_rad;
  
/* Compute new position of punch centre (from velocity * time)
计算冲头中心的新位置（根据速度*时间）*/
   y1=y0 + punch_y_vel * CURRENT_TIME;

/* Set/activate the deforming flag on adjacent cell zone, which      */
/* means that the cells adjacent to the deforming wall will also be  */
/* deformed, in order to avoid skewness.  设置/激活相邻单元格区域上的变形标记，
这意味着与变形壁相邻的单元格也将变形，以避免偏斜。                           */
  SET_DEFORMING_THREAD_FLAG (THREAD_T0 (tf));

/* Loop over the deforming boundary zone's faces; 在变形边界区域的面上循环；                   */
/* inner loop loops over all nodes of a given face;   内部循环遍历给定面的所有节点；               */
/* Thus, since one node can belong to several faces, one must guard 因此，由于一个节点可以属于多个面，因此一个节点必须保护 */
/* against operating on a given node more than once:     禁止在给定节点上多次操作            */

  begin_f_loop (f, tf)
    {
      f_node_loop (f, tf, n)
        {
          node_p = F_NODE (f, tf, n);

      /*We have now identified a given node on the face, but lets check a few things  我们现在已经确定了面上的给定节点，
      但是让我们检查一下几件事 */
      /*Check 1:  Has this node previously been updated?   检查1：此节点以前是否已更新？ */ 
          if (NODE_POS_NEED_UPDATE (node_p))
            {
               NODE_POS_UPDATED (node_p);

      /*Check 2:  Does this node lie within the diameter of the punch?检查2：此节点是否位于冲头的直径之内？ */
               x= NODE_X (node_p);
               if (fabs(x0-x) < circ_rad)
                  {
      /* Compute y-coordinate of punch at this x co-ordinate 在此x坐标处计算冲头的y坐标 */
                   y2= y1 - sqrt( (circ_rad*circ_rad) - (x0-x)*(x0-x) );
      /*Check 3:  If the punch is moving the top plate here, then change node position  检查3：如果打孔机将顶板移到此处，则更改节点位置 */ 
                   if (y2 < 1.0 )
                      {
                         NODE_Y (node_p) = y2;
                      }    /*End check 3 */
                   }       /*End check 2 */
             }             /*End check 1 */

        }       /*End Node Loop*/
    }    /*End Face Loop*/
  end_f_loop (f, tf)
}
/**End of UDF**/
