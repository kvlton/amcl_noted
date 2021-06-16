/*
 *  Player - One Hell of a Robot Server
 *  Copyright (C) 2000  Brian Gerkey   &  Kasper Stoy
 *                      gerkey@usc.edu    kaspers@robotics.usc.edu
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
/**************************************************************************
 * Desc: Range routines
 * Author: Andrew Howard
 * Date: 18 Jan 2003
 * CVS: $Id: map_range.c 1347 2003-05-05 06:24:33Z inspectorg $
**************************************************************************/

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "amcl/map/map.h"

/// 根据单个激光束数据（起点+方向），计算在地图上的理论range（真实距离）
// Extract a single range reading from the map.  Unknown cells and/or
// out-of-bound cells are treated as occupied, which makes it easy to
// use Stage bitmap files.
double map_calc_range(map_t *map, double ox, double oy, double oa, double max_range)
{
  /// Bresenham raytracing 算法
    
  int x0,x1,y0,y1;
  int x,y;
  int xstep, ystep;
  char steep;
  int tmp;
  int deltax, deltay, error, deltaerr;

  // 波束的起点和终点
  x0 = MAP_GXWX(map,ox);
  y0 = MAP_GYWY(map,oy);
  
  x1 = MAP_GXWX(map,ox + max_range * cos(oa));
  y1 = MAP_GYWY(map,oy + max_range * sin(oa));

  // 若y轴差异大于x轴差异，则steep=1
  if(abs(y1-y0) > abs(x1-x0))
    steep = 1;
  else
    steep = 0;

  // 若steep=1，则交换x和y的坐标
  if(steep)
  {
    tmp = x0;
    x0 = y0;
    y0 = tmp;

    tmp = x1;
    x1 = y1;
    y1 = tmp;
  }

  // 计算步长和方向
  deltax = abs(x1-x0);
  deltay = abs(y1-y0);
  error = 0;
  deltaerr = deltay;

  x = x0;
  y = y0;

  if(x0 < x1)
    xstep = 1;
  else
    xstep = -1;
  if(y0 < y1)
    ystep = 1;
  else
    ystep = -1;

  
  // 检查机器人的位置坐标，确保它在地图的上，并且没有处于未知区域或者占用删格中
  if(steep)
  {
    if(!MAP_VALID(map,y,x) || map->cells[MAP_INDEX(map,y,x)].occ_state > -1)
      return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * map->scale;
  }
  else
  {
    if(!MAP_VALID(map,x,y) || map->cells[MAP_INDEX(map,x,y)].occ_state > -1)
      return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * map->scale;
  }

  // 沿着波束延伸的方向依次检查所覆盖的删格是否空闲
  // 如果扫描到波束的终点都是空闲的，我们将返回传感器的最大量程
  // 否则，扫描波束要么触碰到了地图的边界，要么被未知区域或者占用删格阻挡了去路
  while(x != (x1 + xstep * 1))
  {
    // 当前格子的xy坐标
    x += xstep;
    error += deltaerr;
    if(2*error >= deltax)
    {
      y += ystep;
      error -= deltax;
    }

    // 是否是 地图边界、未知区域 或者 占用栅格
    if(steep)
    {
      if(!MAP_VALID(map,y,x) || map->cells[MAP_INDEX(map,y,x)].occ_state > -1)
        return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * map->scale; // 返回最大距离
    }
    else
    {
      if(!MAP_VALID(map,x,y) || map->cells[MAP_INDEX(map,x,y)].occ_state > -1)
        return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * map->scale;
    }
  }
  
  return max_range; // 扫描到波束的终点都是空闲的，返回传感器的最大量程
}
