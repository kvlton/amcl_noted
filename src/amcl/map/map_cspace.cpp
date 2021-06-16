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

#include <queue>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "amcl/map/map.h"

class CellData
{
  public:
    map_t* map_;
    unsigned int i_, j_;         // 当前坐标
    unsigned int src_i_, src_j_; // 起点坐标
};

/// DistanceMap:边长为 max_dist 的正方形
/// 每个格子保存着到原点的距离，通过查表的方式来避免乘法运算
class CachedDistanceMap
{
  public:
    
    /// 构造函数
    CachedDistanceMap(double scale, double max_dist) : 
      distances_(NULL), scale_(scale), max_dist_(max_dist) 
    {
      cell_radius_ = max_dist / scale;
      
      // 边长实际为 cell_radius_+2
      distances_ = new double *[cell_radius_+2];
      for(int i=0; i<=cell_radius_+1; i++)
      {
        distances_[i] = new double[cell_radius_+2];
        for(int j=0; j<=cell_radius_+1; j++)
        {
          distances_[i][j] = sqrt(i*i + j*j); // 到原点的距离
        }
      }
    }
    
    ~CachedDistanceMap()
    {
      if(distances_)
      {
        for(int i=0; i<=cell_radius_+1; i++)
          delete[] distances_[i];
        delete[] distances_;
      }
    }
    
    double** distances_; // 到原点的距离
    double scale_;
    double max_dist_; // 地图的最大扩张距离（m）
    int cell_radius_; // 地图的最大扩张距离（cell）
};


bool operator<(const CellData& a, const CellData& b)
{
  return a.map_->cells[MAP_INDEX(a.map_, a.i_, a.j_)].occ_dist > a.map_->cells[MAP_INDEX(b.map_, b.i_, b.j_)].occ_dist;
}

/// 生成一张 DistanceMap，通过查表代替乘法运算
CachedDistanceMap*
get_distance_map(double scale, double max_dist)
{
  static CachedDistanceMap* cdm = NULL;

  if(!cdm || (cdm->scale_ != scale) || (cdm->max_dist_ != max_dist))
  {
    if(cdm)
      delete cdm;
    cdm = new CachedDistanceMap(scale, max_dist);
  }

  return cdm;
}

/// 加入队列
/// 先更新dist，再加入队列（优先队列会使得队首的dist最小）
void enqueue(map_t* map, int i, int j,
	     int src_i, int src_j,
	     std::priority_queue<CellData>& Q,
	     CachedDistanceMap* cdm,
	     unsigned char* marked)
{
  // 已经加入过队列中
  if(marked[MAP_INDEX(map, i, j)])
    return;

  // 到起点的距离
  int di = abs(i - src_i);
  int dj = abs(j - src_j);
  double distance = cdm->distances_[di][dj]; // 通过查表来代替乘法

  if(distance > cdm->cell_radius_)
    return;

  // 直接更新dist
  map->cells[MAP_INDEX(map, i, j)].occ_dist = distance * map->scale;

  // 加入队列中
  CellData cell;
  cell.map_ = map;
  cell.i_ = i;
  cell.j_ = j;
  cell.src_i_ = src_i;
  cell.src_j_ = src_j;

  Q.push(cell);

  marked[MAP_INDEX(map, i, j)] = 1;
}

/// 似然场预计算（只做一次）
// Update the cspace distance values
void map_update_cspace(map_t *map, double max_occ_dist)
{
  unsigned char* marked; // 是否加入过队列中
  std::priority_queue<CellData> Q;

  marked = new unsigned char[map->size_x*map->size_y];
  memset(marked, 0, sizeof(unsigned char) * map->size_x*map->size_y);

  map->max_occ_dist = max_occ_dist;

  // 生成一张 DistanceMap，通过查表代替乘法运算
  CachedDistanceMap* cdm = get_distance_map(map->scale, map->max_occ_dist);

  // 遍历地图的格子，将障碍物加入到队列中
  // Enqueue all the obstacle cells
  CellData cell;
  cell.map_ = map;
  for(int i=0; i<map->size_x; i++)
  {
    cell.src_i_ = cell.i_ = i; // 障碍物为膨胀的起点
    for(int j=0; j<map->size_y; j++)
    {
      if(map->cells[MAP_INDEX(map, i, j)].occ_state == +1)
      {
        // 将障碍物加入到队列中
        map->cells[MAP_INDEX(map, i, j)].occ_dist = 0.0; // 障碍物的dist预设为0
        cell.src_j_ = cell.j_ = j;
        marked[MAP_INDEX(map, i, j)] = 1;
        Q.push(cell);
      }
      else
        map->cells[MAP_INDEX(map, i, j)].occ_dist = max_occ_dist; // free的dist预设为最大
    }
  }

  // BFS更新地图的dist
  while(!Q.empty())
  {
    // 取出队列最前面的格子（优先队列会保证该格子的dist是最小的）
    CellData current_cell = Q.top();
    
    // 分别将上下左右四个格子加入到队列中（先更新dist，再加入队列）
    if(current_cell.i_ > 0)
      enqueue(map, current_cell.i_-1, current_cell.j_, 
	      current_cell.src_i_, current_cell.src_j_,
	      Q, cdm, marked);
    if(current_cell.j_ > 0)
      enqueue(map, current_cell.i_, current_cell.j_-1, 
	      current_cell.src_i_, current_cell.src_j_,
	      Q, cdm, marked);
    if((int)current_cell.i_ < map->size_x - 1)
      enqueue(map, current_cell.i_+1, current_cell.j_, 
	      current_cell.src_i_, current_cell.src_j_,
	      Q, cdm, marked);
    if((int)current_cell.j_ < map->size_y - 1)
      enqueue(map, current_cell.i_, current_cell.j_+1, 
	      current_cell.src_i_, current_cell.src_j_,
	      Q, cdm, marked);

    Q.pop();
  }

  delete[] marked;
}

#if 0
// TODO: replace this with a more efficient implementation.  Not crucial,
// because we only do it once, at startup.
void map_update_cspace(map_t *map, double max_occ_dist)
{
  int i, j;
  int ni, nj;
  int s;
  double d;
  map_cell_t *cell, *ncell;

  map->max_occ_dist = max_occ_dist;
  s = (int) ceil(map->max_occ_dist / map->scale);

  // Reset the distance values
  for (j = 0; j < map->size_y; j++)
  {
    for (i = 0; i < map->size_x; i++)
    {
      cell = map->cells + MAP_INDEX(map, i, j);
      cell->occ_dist = map->max_occ_dist;
    }
  }

  // Find all the occupied cells and update their neighbours
  for (j = 0; j < map->size_y; j++)
  {
    for (i = 0; i < map->size_x; i++)
    {
      cell = map->cells + MAP_INDEX(map, i, j);
      if (cell->occ_state != +1)
        continue;
          
      cell->occ_dist = 0;

      // Update adjacent cells
      for (nj = -s; nj <= +s; nj++)
      {
        for (ni = -s; ni <= +s; ni++)
        {
          if (!MAP_VALID(map, i + ni, j + nj))
            continue;

          ncell = map->cells + MAP_INDEX(map, i + ni, j + nj);
          d = map->scale * sqrt(ni * ni + nj * nj);

          if (d < ncell->occ_dist)
            ncell->occ_dist = d;
        }
      }
    }
  }
  
  return;
}
#endif
