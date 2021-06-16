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
 * Desc: Simple particle filter for localization.
 * Author: Andrew Howard
 * Date: 10 Dec 2002
 * CVS: $Id: pf.h 3293 2005-11-19 08:37:45Z gerkey $
 *************************************************************************/

#ifndef PF_H
#define PF_H

#include "pf_vector.h"
#include "pf_kdtree.h"

#ifdef __cplusplus
extern "C" {
#endif

// Forward declarations
struct _pf_t;
struct _rtk_fig_t;
struct _pf_sample_set_t;

/// 回调函数：随机生成粒子
// Function prototype for the initialization model; generates a sample pose from
// an appropriate distribution.
typedef pf_vector_t (*pf_init_model_fn_t) (void *init_data);

/// 回调函数：运动模型（unused）
// Function prototype for the action model; generates a sample pose from
// an appropriate distribution
typedef void (*pf_action_model_fn_t) (void *action_data, 
                                      struct _pf_sample_set_t* set);

/// 回调函数：测量模型，根据激光数据计算粒子的权重
// Function prototype for the sensor model; determines the probability
// for the given set of sample poses.
typedef double (*pf_sensor_model_fn_t) (void *sensor_data, 
                                        struct _pf_sample_set_t* set);


/// 单个粒子
// Information for a single sample
typedef struct
{
  // Pose represented by this sample
  pf_vector_t pose;

  // Weight for this pose
  double weight;
  
} pf_sample_t;


/// 粒子簇，用于描述多峰分布，一个粒子簇表示一个峰
// Information for a cluster of samples
typedef struct
{
  // Number of samples
  int count;     // 粒子数量

  // Total weight of samples in this cluster
  double weight; // 粒子总权重

  // Cluster statistics
  pf_vector_t mean;
  pf_matrix_t cov;

  // Workspace
  double m[4], c[2][2]; // 均值（x y cosθ sinθ）y、协方差（x y），用于计算 mean cov
  
} pf_cluster_t;


/// 粒子集合
// Information for a set of samples
typedef struct _pf_sample_set_t
{
  // The samples
  int sample_count; // 粒子数
  pf_sample_t *samples;

  // A kdtree encoding the histogram
  pf_kdtree_t *kdtree;

  // Clusters
  int cluster_count, cluster_max_count;
  pf_cluster_t *clusters;

  // Filter statistics
  pf_vector_t mean;
  pf_matrix_t cov;
  int converged; 
  double n_effective; // 粒子权重的集中程度
} pf_sample_set_t;


/// 滤波器
// Information for an entire filter
typedef struct _pf_t
{
  // This min and max number of samples
  int min_samples, max_samples; // 重采样的最小粒子数和最大粒子数

  // Population size parameters
  double pop_err, pop_z; // KLD采样的两个参数

  // Resample limit cache
  // 当直方图的bin数k确定时，对应的重采样粒子数Mx也确定
  int *limit_cache; // 保存 KLD算法 计算得到的 重采样粒子数 Mx，避免重复计算
  
  // The sample sets.  We keep two sets and use [current_set]
  // to identify the active set.
  int current_set; // 当前粒子集的索引 ＡctiveSet
  pf_sample_set_t sets[2]; // 一个作为ＡctiveＳet 一个用来保存采样后的粒子（重采样后的ＡctiveSet）

  // Running averages, slow and fast, of likelihood
  double w_slow, w_fast; // 用于随机生成粒子

  // Decay rates for running averages
  double alpha_slow, alpha_fast;

  // Function used to draw random pose samples
  pf_init_model_fn_t random_pose_fn; // 回调函数 AmclNode::uniformPoseGenerator 生成随机粒子
  void *random_pose_data; // 实际上是 map_t 数据结构

  double dist_threshold; //distance threshold in each axis over which the pf is considered to not be converged
  int converged; 

  // boolean parameter to enamble/diable selective resampling
  int selective_resampling;  // 有选择性地重采样（不一定要进行重采样）
} pf_t;

/// 创建滤波器
// Create a new filter
pf_t *pf_alloc(int min_samples, int max_samples,
               double alpha_slow, double alpha_fast,
               pf_init_model_fn_t random_pose_fn, void *random_pose_data);

/// 释放滤波器
// Free an existing filter
void pf_free(pf_t *pf);

/// 用高斯分布来初始化滤波器
// Initialize the filter using a guassian
void pf_init(pf_t *pf, pf_vector_t mean, pf_matrix_t cov);

/// 初始化滤波器（使用 自定义模型 代替 高斯模型）
// Initialize the filter using some model
void pf_init_model(pf_t *pf, pf_init_model_fn_t init_fn, void *init_data);

/// 运动更新 未使用到（在odom中进行运动更新）
// Update the filter with some new action
void pf_update_action(pf_t *pf, pf_action_model_fn_t action_fn, void *action_data);

/// 测量更新 使用激光传感器数据来更新粒子滤波器
// Update the filter with some new sensor observation
void pf_update_sensor(pf_t *pf, pf_sensor_model_fn_t sensor_fn, void *sensor_data);

/// 重采样
// Resample the distribution
void pf_update_resample(pf_t *pf);

// set selective resampling parameter
void pf_set_selective_resampling(pf_t *pf, int selective_resampling);

/// 计算CEP统计，for pf_draw
// Compute the CEP statistics (mean and variance).
void pf_get_cep_stats(pf_t *pf, pf_vector_t *mean, double *var);

/// 获取某一聚类的统计特性
// Compute the statistics for a particular cluster.  Returns 0 if
// there is no such cluster.
int pf_get_cluster_stats(pf_t *pf, int cluster, double *weight,
                         pf_vector_t *mean, pf_matrix_t *cov);

/// 更新粒子集的粒子簇
// Re-compute the cluster statistics for a sample set
void pf_cluster_stats(pf_t *pf, pf_sample_set_t *set);


// Display the sample set
void pf_draw_samples(pf_t *pf, struct _rtk_fig_t *fig, int max_samples);

// Draw the histogram (kdtree)
void pf_draw_hist(pf_t *pf, struct _rtk_fig_t *fig);

// Draw the CEP statistics
void pf_draw_cep_stats(pf_t *pf, struct _rtk_fig_t *fig);

// Draw the cluster statistics
void pf_draw_cluster_stats(pf_t *pf, struct _rtk_fig_t *fig);

/// 检查滤波器是否收敛 
//calculate if the particle filter has converged - 
//and sets the converged flag in the current set and the pf 
int pf_update_converged(pf_t *pf);

/// 初始化收敛状态 
//sets the current set and pf converged values to zero
void pf_init_converged(pf_t *pf);

void pf_copy_set(pf_sample_set_t* set_a, pf_sample_set_t* set_b);

#ifdef __cplusplus
}
#endif


#endif
