-- MySQL dump 10.11
--
-- Host: mcs10    Database: vrtrack_trash
-- ------------------------------------------------------
-- Server version	5.0.37-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `allocation`
--

DROP TABLE IF EXISTS `allocation`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `allocation` (
  `study_id` smallint(5) unsigned NOT NULL default '0',
  `individual_id` smallint(5) unsigned NOT NULL default '0',
  `seq_centre_id` smallint(5) unsigned NOT NULL default '0',
  PRIMARY KEY  (`study_id`,`individual_id`,`seq_centre_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `allocation`
--

LOCK TABLES `allocation` WRITE;
/*!40000 ALTER TABLE `allocation` DISABLE KEYS */;
/*!40000 ALTER TABLE `allocation` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `assembly`
--

DROP TABLE IF EXISTS `assembly`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `assembly` (
  `assembly_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `reference_size` int(11) default NULL,
  `taxon_id` mediumint(8) unsigned default NULL,
  `translation_table` smallint(5) unsigned default NULL,
  PRIMARY KEY  (`assembly_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `assembly`
--

LOCK TABLES `assembly` WRITE;
/*!40000 ALTER TABLE `assembly` DISABLE KEYS */;
/*!40000 ALTER TABLE `assembly` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `autoqc`
--

DROP TABLE IF EXISTS `autoqc`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `autoqc` (
  `autoqc_id` mediumint(8) unsigned NOT NULL auto_increment,
  `mapstats_id` mediumint(8) unsigned NOT NULL default '0',
  `test` varchar(50) NOT NULL default '',
  `result` tinyint(1) default '0',
  `reason` varchar(200) NOT NULL default '',
  PRIMARY KEY  (`autoqc_id`),
  UNIQUE KEY `mapstats_test` (`mapstats_id`,`test`),
  KEY `mapstats_id` (`mapstats_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `autoqc`
--

LOCK TABLES `autoqc` WRITE;
/*!40000 ALTER TABLE `autoqc` DISABLE KEYS */;
/*!40000 ALTER TABLE `autoqc` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `exome_design`
--

DROP TABLE IF EXISTS `exome_design`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `exome_design` (
  `exome_design_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `bait_bases` bigint(20) unsigned default NULL,
  `target_bases` bigint(20) unsigned default NULL,
  PRIMARY KEY  (`exome_design_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `exome_design`
--

LOCK TABLES `exome_design` WRITE;
/*!40000 ALTER TABLE `exome_design` DISABLE KEYS */;
/*!40000 ALTER TABLE `exome_design` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `file`
--

DROP TABLE IF EXISTS `file`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `file` (
  `row_id` int(10) unsigned NOT NULL auto_increment,
  `file_id` mediumint(8) unsigned NOT NULL default '0',
  `lane_id` mediumint(8) unsigned NOT NULL default '0',
  `name` varchar(255) NOT NULL default '',
  `hierarchy_name` varchar(255) default NULL,
  `processed` int(10) default '0',
  `type` tinyint(4) default NULL,
  `readlen` smallint(5) unsigned default NULL,
  `raw_reads` bigint(20) unsigned default NULL,
  `raw_bases` bigint(20) unsigned default NULL,
  `mean_q` float unsigned default NULL,
  `md5` char(32) default NULL,
  `note_id` mediumint(8) unsigned default NULL,
  `changed` datetime NOT NULL default '0000-00-00 00:00:00',
  `latest` tinyint(1) default '0',
  PRIMARY KEY  (`row_id`),
  KEY `file_id` (`file_id`),
  KEY `lane_id` (`lane_id`),
  KEY `hierarchy_name` (`hierarchy_name`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `file`
--

LOCK TABLES `file` WRITE;
/*!40000 ALTER TABLE `file` DISABLE KEYS */;
INSERT INTO `file` VALUES (98,97,97,'9252616016_K_Grn.idat','9252616016_K_Grn.idat',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-08-09 14:07:57',0),(99,97,97,'9252616016_K_Grn.idat','9252616016_K_Grn.idat',0,8,NULL,NULL,NULL,NULL,'45fe61ebdb49e343558e36aff8dd20b9',NULL,'2013-08-09 14:07:57',1),(101,100,100,'9252616016_I_Grn.idat','9252616016_I_Grn.idat',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-08-09 14:07:57',0),(102,100,100,'9252616016_I_Grn.idat','9252616016_I_Grn.idat',0,8,NULL,NULL,NULL,NULL,'ae61b62dc73776dd8af04691655e9037',NULL,'2013-08-09 14:07:57',1),(104,103,103,'9252616016_G_Grn.idat','9252616016_G_Grn.idat',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-08-09 14:07:57',0),(105,103,103,'9252616016_G_Grn.idat','9252616016_G_Grn.idat',0,8,NULL,NULL,NULL,NULL,'e2f7a2635d5225d36914648f0ddfad9f',NULL,'2013-08-09 14:07:57',1),(113,112,112,'9252616016_B_Grn.idat','9252616016_B_Grn.idat',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-08-09 14:07:58',0),(114,112,112,'9252616016_B_Grn.idat','9252616016_B_Grn.idat',0,8,NULL,NULL,NULL,NULL,'f9ee6f9238449ef05caebdfc01b52e4a',NULL,'2013-08-09 14:07:58',1);
/*!40000 ALTER TABLE `file` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `image`
--

DROP TABLE IF EXISTS `image`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `image` (
  `image_id` mediumint(8) unsigned NOT NULL auto_increment,
  `mapstats_id` mediumint(8) unsigned NOT NULL default '0',
  `name` varchar(255) NOT NULL default '',
  `caption` varchar(40) default NULL,
  `image` mediumblob,
  PRIMARY KEY  (`image_id`),
  KEY `mapstats_id` (`mapstats_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `image`
--

LOCK TABLES `image` WRITE;
/*!40000 ALTER TABLE `image` DISABLE KEYS */;
/*!40000 ALTER TABLE `image` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `individual`
--

DROP TABLE IF EXISTS `individual`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `individual` (
  `individual_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `hierarchy_name` varchar(255) NOT NULL default '',
  `alias` varchar(40) NOT NULL default '',
  `sex` enum('M','F','unknown') default 'unknown',
  `acc` varchar(40) default NULL,
  `species_id` smallint(5) unsigned default NULL,
  `population_id` smallint(5) unsigned default NULL,
  PRIMARY KEY  (`individual_id`),
  UNIQUE KEY `name` (`name`),
  UNIQUE KEY `hierarchy_name` (`hierarchy_name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `individual`
--

LOCK TABLES `individual` WRITE;
/*!40000 ALTER TABLE `individual` DISABLE KEYS */;
INSERT INTO `individual` VALUES (9,'27af9a9b-01b2-4cb6-acef-ea52d83e3d26','27af9a9b_01b2_4cb6_acef_ea52d83e3d26','','unknown',NULL,1,1);
/*!40000 ALTER TABLE `individual` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `lane`
--

DROP TABLE IF EXISTS `lane`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `lane` (
  `row_id` int(10) unsigned NOT NULL auto_increment,
  `lane_id` mediumint(8) unsigned NOT NULL default '0',
  `library_id` smallint(5) unsigned NOT NULL default '0',
  `seq_request_id` mediumint(8) unsigned NOT NULL default '0',
  `name` varchar(255) NOT NULL default '',
  `hierarchy_name` varchar(255) NOT NULL default '',
  `acc` varchar(40) default NULL,
  `readlen` smallint(5) unsigned default NULL,
  `paired` tinyint(1) default NULL,
  `raw_reads` bigint(20) unsigned default NULL,
  `raw_bases` bigint(20) unsigned default NULL,
  `npg_qc_status` enum('pending','pass','fail','-') default 'pending',
  `processed` int(10) default '0',
  `auto_qc_status` enum('no_qc','passed','failed') default 'no_qc',
  `qc_status` enum('no_qc','pending','passed','failed','gt_pending','investigate') default 'no_qc',
  `gt_status` enum('unchecked','confirmed','wrong','unconfirmed','candidate','unknown','swapped') default 'unchecked',
  `submission_id` smallint(5) unsigned default NULL,
  `withdrawn` tinyint(1) default NULL,
  `note_id` mediumint(8) unsigned default NULL,
  `changed` datetime NOT NULL default '0000-00-00 00:00:00',
  `run_date` datetime default NULL,
  `storage_path` varchar(255) default NULL,
  `latest` tinyint(1) default '0',
  PRIMARY KEY  (`row_id`),
  KEY `lane_id` (`lane_id`),
  KEY `lanename` (`name`),
  KEY `library_id` (`library_id`),
  KEY `hierarchy_name` (`hierarchy_name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `lane`
--

LOCK TABLES `lane` WRITE;
/*!40000 ALTER TABLE `lane` DISABLE KEYS */;
INSERT INTO `lane` VALUES (97,97,0,0,'9252616016_K_Grn','9252616016_K_Grn',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,NULL,0),(98,97,129,0,'9252616016_K_Grn','9252616016_K_Grn',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,NULL,0),(99,97,129,0,'9252616016_K_Grn','9252616016_K_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',0),(100,100,0,0,'9252616016_I_Grn','9252616016_I_Grn',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,NULL,0),(101,100,133,0,'9252616016_I_Grn','9252616016_I_Grn',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,NULL,0),(102,100,133,0,'9252616016_I_Grn','9252616016_I_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',0),(103,103,0,0,'9252616016_G_Grn','9252616016_G_Grn',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,NULL,0),(104,103,137,0,'9252616016_G_Grn','9252616016_G_Grn',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,NULL,0),(105,103,137,0,'9252616016_G_Grn','9252616016_G_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:57',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',0),(112,112,0,0,'9252616016_B_Grn','9252616016_B_Grn',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:58',NULL,NULL,0),(113,112,149,0,'9252616016_B_Grn','9252616016_B_Grn',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:58',NULL,NULL,0),(114,112,149,0,'9252616016_B_Grn','9252616016_B_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-09 14:07:58',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',0),(144,97,129,0,'9252616016_K_Grn','9252616016_K_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',1,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-15 16:54:04',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',0),(145,97,129,0,'9252616016_K_Grn','9252616016_K_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',5,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-15 16:54:04',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',1),(150,112,149,0,'9252616016_B_Grn','9252616016_B_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',1,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-15 16:54:04',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',0),(151,112,149,0,'9252616016_B_Grn','9252616016_B_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',5,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-15 16:54:04',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',1),(178,103,137,0,'9252616016_G_Grn','9252616016_G_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',1,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-15 16:54:11',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',0),(179,103,137,0,'9252616016_G_Grn','9252616016_G_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',5,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-15 16:54:11',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',1),(186,100,133,0,'9252616016_I_Grn','9252616016_I_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',1,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-15 16:54:16',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',0),(187,100,133,0,'9252616016_I_Grn','9252616016_I_Grn','cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',NULL,0,0,0,'pending',5,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-08-15 16:54:16',NULL,'/lustre/scratch105/vrpipe/refs/human/ncbi37/resources_hipsci/expression/downloads/cf7095aa-363e-43aa-8c85-09fdc6ffc9cb',1);
/*!40000 ALTER TABLE `lane` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Temporary table structure for view `latest_file`
--

DROP TABLE IF EXISTS `latest_file`;
/*!50001 DROP VIEW IF EXISTS `latest_file`*/;
/*!50001 CREATE TABLE `latest_file` (
  `row_id` int(10) unsigned,
  `file_id` mediumint(8) unsigned,
  `lane_id` mediumint(8) unsigned,
  `name` varchar(255),
  `hierarchy_name` varchar(255),
  `processed` int(10),
  `type` tinyint(4),
  `readlen` smallint(5) unsigned,
  `raw_reads` bigint(20) unsigned,
  `raw_bases` bigint(20) unsigned,
  `mean_q` float unsigned,
  `md5` char(32),
  `note_id` mediumint(8) unsigned,
  `changed` datetime,
  `latest` tinyint(1)
) */;

--
-- Temporary table structure for view `latest_lane`
--

DROP TABLE IF EXISTS `latest_lane`;
/*!50001 DROP VIEW IF EXISTS `latest_lane`*/;
/*!50001 CREATE TABLE `latest_lane` (
  `row_id` int(10) unsigned,
  `lane_id` mediumint(8) unsigned,
  `library_id` smallint(5) unsigned,
  `seq_request_id` mediumint(8) unsigned,
  `name` varchar(255),
  `hierarchy_name` varchar(255),
  `acc` varchar(40),
  `readlen` smallint(5) unsigned,
  `paired` tinyint(1),
  `raw_reads` bigint(20) unsigned,
  `raw_bases` bigint(20) unsigned,
  `npg_qc_status` enum('pending','pass','fail','-'),
  `processed` int(10),
  `auto_qc_status` enum('no_qc','passed','failed'),
  `qc_status` enum('no_qc','pending','passed','failed','gt_pending','investigate'),
  `gt_status` enum('unchecked','confirmed','wrong','unconfirmed','candidate','unknown','swapped'),
  `submission_id` smallint(5) unsigned,
  `withdrawn` tinyint(1),
  `note_id` mediumint(8) unsigned,
  `changed` datetime,
  `run_date` datetime,
  `storage_path` varchar(255),
  `latest` tinyint(1)
) */;

--
-- Temporary table structure for view `latest_library`
--

DROP TABLE IF EXISTS `latest_library`;
/*!50001 DROP VIEW IF EXISTS `latest_library`*/;
/*!50001 CREATE TABLE `latest_library` (
  `row_id` int(10) unsigned,
  `library_id` smallint(5) unsigned,
  `library_request_id` mediumint(8) unsigned,
  `sample_id` smallint(5) unsigned,
  `ssid` mediumint(8) unsigned,
  `name` varchar(255),
  `hierarchy_name` varchar(255),
  `prep_status` enum('unknown','pending','started','passed','failed','cancelled','hold'),
  `auto_qc_status` enum('no_qc','passed','failed'),
  `qc_status` enum('no_qc','pending','passed','failed'),
  `fragment_size_from` mediumint(8) unsigned,
  `fragment_size_to` mediumint(8) unsigned,
  `library_type_id` smallint(5) unsigned,
  `library_tag` smallint(5) unsigned,
  `library_tag_group` smallint(5) unsigned,
  `library_tag_sequence` varchar(1024),
  `seq_centre_id` smallint(5) unsigned,
  `seq_tech_id` smallint(5) unsigned,
  `open` tinyint(1),
  `note_id` mediumint(8) unsigned,
  `changed` datetime,
  `latest` tinyint(1)
) */;

--
-- Temporary table structure for view `latest_library_request`
--

DROP TABLE IF EXISTS `latest_library_request`;
/*!50001 DROP VIEW IF EXISTS `latest_library_request`*/;
/*!50001 CREATE TABLE `latest_library_request` (
  `row_id` int(10) unsigned,
  `library_request_id` mediumint(8) unsigned,
  `sample_id` smallint(5) unsigned,
  `ssid` mediumint(8) unsigned,
  `prep_status` enum('unknown','pending','started','passed','failed','cancelled','hold'),
  `note_id` mediumint(8) unsigned,
  `changed` datetime,
  `latest` tinyint(1)
) */;

--
-- Temporary table structure for view `latest_mapstats`
--

DROP TABLE IF EXISTS `latest_mapstats`;
/*!50001 DROP VIEW IF EXISTS `latest_mapstats`*/;
/*!50001 CREATE TABLE `latest_mapstats` (
  `row_id` int(10) unsigned,
  `mapstats_id` mediumint(8) unsigned,
  `lane_id` mediumint(8) unsigned,
  `mapper_id` smallint(5) unsigned,
  `assembly_id` smallint(5) unsigned,
  `raw_reads` bigint(20) unsigned,
  `raw_bases` bigint(20) unsigned,
  `clip_bases` bigint(20) unsigned,
  `reads_mapped` bigint(20) unsigned,
  `reads_paired` bigint(20) unsigned,
  `bases_mapped` bigint(20) unsigned,
  `rmdup_reads_mapped` bigint(20) unsigned,
  `rmdup_bases_mapped` bigint(20) unsigned,
  `adapter_reads` bigint(20) unsigned,
  `error_rate` float unsigned,
  `mean_insert` float unsigned,
  `sd_insert` float unsigned,
  `gt_expected` varchar(40),
  `gt_found` varchar(40),
  `gt_ratio` float unsigned,
  `note_id` mediumint(8) unsigned,
  `changed` datetime,
  `latest` tinyint(1),
  `bait_near_bases_mapped` bigint(20) unsigned,
  `target_near_bases_mapped` bigint(20) unsigned,
  `bait_bases_mapped` bigint(20) unsigned,
  `mean_bait_coverage` float unsigned,
  `bait_coverage_sd` float unsigned,
  `off_bait_bases` bigint(20) unsigned,
  `reads_on_bait` bigint(20) unsigned,
  `reads_on_bait_near` bigint(20) unsigned,
  `reads_on_target` bigint(20) unsigned,
  `reads_on_target_near` bigint(20) unsigned,
  `target_bases_mapped` bigint(20) unsigned,
  `mean_target_coverage` float unsigned,
  `target_coverage_sd` float unsigned,
  `target_bases_1X` float unsigned,
  `target_bases_2X` float unsigned,
  `target_bases_5X` float unsigned,
  `target_bases_10X` float unsigned,
  `target_bases_20X` float unsigned,
  `target_bases_50X` float unsigned,
  `target_bases_100X` float unsigned,
  `exome_design_id` smallint(5) unsigned,
  `percentage_reads_with_transposon` float unsigned,
  `is_qc` tinyint(1),
  `prefix` varchar(40)
) */;

--
-- Temporary table structure for view `latest_project`
--

DROP TABLE IF EXISTS `latest_project`;
/*!50001 DROP VIEW IF EXISTS `latest_project`*/;
/*!50001 CREATE TABLE `latest_project` (
  `row_id` int(10) unsigned,
  `project_id` smallint(5) unsigned,
  `ssid` mediumint(8) unsigned,
  `name` varchar(255),
  `hierarchy_name` varchar(255),
  `study_id` smallint(5),
  `note_id` mediumint(8) unsigned,
  `changed` datetime,
  `latest` tinyint(1)
) */;

--
-- Temporary table structure for view `latest_sample`
--

DROP TABLE IF EXISTS `latest_sample`;
/*!50001 DROP VIEW IF EXISTS `latest_sample`*/;
/*!50001 CREATE TABLE `latest_sample` (
  `row_id` int(10) unsigned,
  `sample_id` smallint(5) unsigned,
  `project_id` smallint(5) unsigned,
  `ssid` mediumint(8) unsigned,
  `name` varchar(255),
  `hierarchy_name` varchar(40),
  `individual_id` smallint(5) unsigned,
  `note_id` mediumint(8) unsigned,
  `changed` datetime,
  `latest` tinyint(1)
) */;

--
-- Temporary table structure for view `latest_seq_request`
--

DROP TABLE IF EXISTS `latest_seq_request`;
/*!50001 DROP VIEW IF EXISTS `latest_seq_request`*/;
/*!50001 CREATE TABLE `latest_seq_request` (
  `row_id` int(10) unsigned,
  `seq_request_id` mediumint(8) unsigned,
  `library_id` smallint(5) unsigned,
  `multiplex_pool_id` smallint(5) unsigned,
  `ssid` mediumint(8) unsigned,
  `seq_type` enum('HiSeq Paired end sequencing','Illumina-A HiSeq Paired end sequencing','Illumina-A Paired end sequencing','Illumina-A Pulldown ISC','Illumina-A Pulldown SC','Illumina-A Pulldown WGS','Illumina-A Single ended hi seq sequencing','Illumina-A Single ended sequencing','Illumina-B HiSeq Paired end sequencing','Illumina-B Paired end sequencing','Illumina-B Single ended hi seq sequencing','Illumina-B Single ended sequencing','Illumina-C HiSeq Paired end sequencing','Illumina-C MiSeq sequencing','Illumina-C Paired end sequencing','Illumina-C Single ended hi seq sequencing','Illumina-C Single ended sequencing','MiSeq sequencing','Paired end sequencing','Single ended hi seq sequencing','Single Ended Hi Seq Sequencing Control','Single ended sequencing'),
  `seq_status` enum('unknown','pending','started','passed','failed','cancelled','hold'),
  `note_id` mediumint(8) unsigned,
  `changed` datetime,
  `latest` tinyint(1)
) */;

--
-- Table structure for table `library`
--

DROP TABLE IF EXISTS `library`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `library` (
  `row_id` int(10) unsigned NOT NULL auto_increment,
  `library_id` smallint(5) unsigned NOT NULL default '0',
  `library_request_id` mediumint(8) unsigned NOT NULL default '0',
  `sample_id` smallint(5) unsigned NOT NULL default '0',
  `ssid` mediumint(8) unsigned default NULL,
  `name` varchar(255) NOT NULL default '',
  `hierarchy_name` varchar(255) NOT NULL default '',
  `prep_status` enum('unknown','pending','started','passed','failed','cancelled','hold') default 'unknown',
  `auto_qc_status` enum('no_qc','passed','failed') default 'no_qc',
  `qc_status` enum('no_qc','pending','passed','failed') default 'no_qc',
  `fragment_size_from` mediumint(8) unsigned default NULL,
  `fragment_size_to` mediumint(8) unsigned default NULL,
  `library_type_id` smallint(5) unsigned default NULL,
  `library_tag` smallint(5) unsigned default NULL,
  `library_tag_group` smallint(5) unsigned default NULL,
  `library_tag_sequence` varchar(1024) default NULL,
  `seq_centre_id` smallint(5) unsigned default NULL,
  `seq_tech_id` smallint(5) unsigned default NULL,
  `open` tinyint(1) default '1',
  `note_id` mediumint(8) unsigned default NULL,
  `changed` datetime NOT NULL default '0000-00-00 00:00:00',
  `latest` tinyint(1) default '0',
  PRIMARY KEY  (`row_id`),
  KEY `ssid` (`ssid`),
  KEY `name` (`name`),
  KEY `hierarchy_name` (`hierarchy_name`),
  KEY `sample_id` (`sample_id`),
  KEY `library_id` (`library_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `library`
--

LOCK TABLES `library` WRITE;
/*!40000 ALTER TABLE `library` DISABLE KEYS */;
INSERT INTO `library` VALUES (129,129,0,0,NULL,'qc1hip5529784_9252616016_K','qc1hip5529784_9252616016_K','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(130,129,0,97,NULL,'qc1hip5529784_9252616016_K','qc1hip5529784_9252616016_K','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(131,129,0,97,6016284,'qc1hip5529784_9252616016_K','qc1hip5529784_9252616016_K','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(132,129,0,97,6016284,'qc1hip5529784_9252616016_K','qc1hip5529784_9252616016_K','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'9252616016_K',1,1,1,NULL,'2013-08-09 14:07:57',1),(133,133,0,0,NULL,'qc1hip5529783_9252616016_I','qc1hip5529783_9252616016_I','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(134,133,0,100,NULL,'qc1hip5529783_9252616016_I','qc1hip5529783_9252616016_I','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(135,133,0,100,6016283,'qc1hip5529783_9252616016_I','qc1hip5529783_9252616016_I','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(136,133,0,100,6016283,'qc1hip5529783_9252616016_I','qc1hip5529783_9252616016_I','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'9252616016_I',1,1,1,NULL,'2013-08-09 14:07:57',1),(137,137,0,0,NULL,'qc1hip5529782_9252616016_G','qc1hip5529782_9252616016_G','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(138,137,0,103,NULL,'qc1hip5529782_9252616016_G','qc1hip5529782_9252616016_G','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(139,137,0,103,6016282,'qc1hip5529782_9252616016_G','qc1hip5529782_9252616016_G','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:57',0),(140,137,0,103,6016282,'qc1hip5529782_9252616016_G','qc1hip5529782_9252616016_G','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'9252616016_G',1,1,1,NULL,'2013-08-09 14:07:57',1),(149,149,0,0,NULL,'qc1hip5529785_9252616016_B','qc1hip5529785_9252616016_B','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:58',0),(150,149,0,112,NULL,'qc1hip5529785_9252616016_B','qc1hip5529785_9252616016_B','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:58',0),(151,149,0,112,6016285,'qc1hip5529785_9252616016_B','qc1hip5529785_9252616016_B','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-08-09 14:07:58',0),(152,149,0,112,6016285,'qc1hip5529785_9252616016_B','qc1hip5529785_9252616016_B','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'9252616016_B',1,1,1,NULL,'2013-08-09 14:07:58',1);
/*!40000 ALTER TABLE `library` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `library_multiplex_pool`
--

DROP TABLE IF EXISTS `library_multiplex_pool`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `library_multiplex_pool` (
  `library_multiplex_pool_id` mediumint(8) unsigned NOT NULL auto_increment,
  `multiplex_pool_id` smallint(5) unsigned NOT NULL default '0',
  `library_id` smallint(5) unsigned NOT NULL default '0',
  PRIMARY KEY  (`library_multiplex_pool_id`),
  KEY `library_multiplex_pool_id` (`library_multiplex_pool_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `library_multiplex_pool`
--

LOCK TABLES `library_multiplex_pool` WRITE;
/*!40000 ALTER TABLE `library_multiplex_pool` DISABLE KEYS */;
/*!40000 ALTER TABLE `library_multiplex_pool` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `library_request`
--

DROP TABLE IF EXISTS `library_request`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `library_request` (
  `row_id` int(10) unsigned NOT NULL auto_increment,
  `library_request_id` mediumint(8) unsigned NOT NULL default '0',
  `sample_id` smallint(5) unsigned NOT NULL default '0',
  `ssid` mediumint(8) unsigned default NULL,
  `prep_status` enum('unknown','pending','started','passed','failed','cancelled','hold') default 'unknown',
  `note_id` mediumint(8) unsigned default NULL,
  `changed` datetime NOT NULL default '0000-00-00 00:00:00',
  `latest` tinyint(1) default '0',
  PRIMARY KEY  (`row_id`),
  KEY `library_request_id` (`library_request_id`),
  KEY `ssid` (`ssid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `library_request`
--

LOCK TABLES `library_request` WRITE;
/*!40000 ALTER TABLE `library_request` DISABLE KEYS */;
/*!40000 ALTER TABLE `library_request` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `library_type`
--

DROP TABLE IF EXISTS `library_type`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `library_type` (
  `library_type_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`library_type_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `library_type`
--

LOCK TABLES `library_type` WRITE;
/*!40000 ALTER TABLE `library_type` DISABLE KEYS */;
/*!40000 ALTER TABLE `library_type` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `mapper`
--

DROP TABLE IF EXISTS `mapper`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `mapper` (
  `mapper_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `version` varchar(40) NOT NULL default '0',
  PRIMARY KEY  (`mapper_id`),
  UNIQUE KEY `name_v` (`name`,`version`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `mapper`
--

LOCK TABLES `mapper` WRITE;
/*!40000 ALTER TABLE `mapper` DISABLE KEYS */;
/*!40000 ALTER TABLE `mapper` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `mapstats`
--

DROP TABLE IF EXISTS `mapstats`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `mapstats` (
  `row_id` int(10) unsigned NOT NULL auto_increment,
  `mapstats_id` mediumint(8) unsigned NOT NULL default '0',
  `lane_id` mediumint(8) unsigned NOT NULL default '0',
  `mapper_id` smallint(5) unsigned default NULL,
  `assembly_id` smallint(5) unsigned default NULL,
  `raw_reads` bigint(20) unsigned default NULL,
  `raw_bases` bigint(20) unsigned default NULL,
  `clip_bases` bigint(20) unsigned default NULL,
  `reads_mapped` bigint(20) unsigned default NULL,
  `reads_paired` bigint(20) unsigned default NULL,
  `bases_mapped` bigint(20) unsigned default NULL,
  `rmdup_reads_mapped` bigint(20) unsigned default NULL,
  `rmdup_bases_mapped` bigint(20) unsigned default NULL,
  `adapter_reads` bigint(20) unsigned default NULL,
  `error_rate` float unsigned default NULL,
  `mean_insert` float unsigned default NULL,
  `sd_insert` float unsigned default NULL,
  `gt_expected` varchar(40) default NULL,
  `gt_found` varchar(40) default NULL,
  `gt_ratio` float unsigned default NULL,
  `note_id` mediumint(8) unsigned default NULL,
  `changed` datetime NOT NULL default '0000-00-00 00:00:00',
  `latest` tinyint(1) default '0',
  `bait_near_bases_mapped` bigint(20) unsigned default NULL,
  `target_near_bases_mapped` bigint(20) unsigned default NULL,
  `bait_bases_mapped` bigint(20) unsigned default NULL,
  `mean_bait_coverage` float unsigned default NULL,
  `bait_coverage_sd` float unsigned default NULL,
  `off_bait_bases` bigint(20) unsigned default NULL,
  `reads_on_bait` bigint(20) unsigned default NULL,
  `reads_on_bait_near` bigint(20) unsigned default NULL,
  `reads_on_target` bigint(20) unsigned default NULL,
  `reads_on_target_near` bigint(20) unsigned default NULL,
  `target_bases_mapped` bigint(20) unsigned default NULL,
  `mean_target_coverage` float unsigned default NULL,
  `target_coverage_sd` float unsigned default NULL,
  `target_bases_1X` float unsigned default NULL,
  `target_bases_2X` float unsigned default NULL,
  `target_bases_5X` float unsigned default NULL,
  `target_bases_10X` float unsigned default NULL,
  `target_bases_20X` float unsigned default NULL,
  `target_bases_50X` float unsigned default NULL,
  `target_bases_100X` float unsigned default NULL,
  `exome_design_id` smallint(5) unsigned default NULL,
  `percentage_reads_with_transposon` float unsigned default NULL,
  `is_qc` tinyint(1) default '0',
  `prefix` varchar(40) default '_',
  PRIMARY KEY  (`row_id`),
  KEY `mapstats_id` (`mapstats_id`),
  KEY `lane_id` (`lane_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `mapstats`
--

LOCK TABLES `mapstats` WRITE;
/*!40000 ALTER TABLE `mapstats` DISABLE KEYS */;
/*!40000 ALTER TABLE `mapstats` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `multiplex_pool`
--

DROP TABLE IF EXISTS `multiplex_pool`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `multiplex_pool` (
  `multiplex_pool_id` mediumint(8) unsigned NOT NULL auto_increment,
  `ssid` mediumint(8) unsigned default NULL,
  `name` varchar(255) NOT NULL default '',
  `note_id` mediumint(8) unsigned default NULL,
  PRIMARY KEY  (`multiplex_pool_id`),
  UNIQUE KEY `ssid` (`ssid`),
  KEY `multiplex_pool_id` (`multiplex_pool_id`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `multiplex_pool`
--

LOCK TABLES `multiplex_pool` WRITE;
/*!40000 ALTER TABLE `multiplex_pool` DISABLE KEYS */;
/*!40000 ALTER TABLE `multiplex_pool` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `note`
--

DROP TABLE IF EXISTS `note`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `note` (
  `note_id` mediumint(8) unsigned NOT NULL auto_increment,
  `note` text,
  PRIMARY KEY  (`note_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `note`
--

LOCK TABLES `note` WRITE;
/*!40000 ALTER TABLE `note` DISABLE KEYS */;
/*!40000 ALTER TABLE `note` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `population`
--

DROP TABLE IF EXISTS `population`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `population` (
  `population_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`population_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `population`
--

LOCK TABLES `population` WRITE;
/*!40000 ALTER TABLE `population` DISABLE KEYS */;
INSERT INTO `population` VALUES (1,'Population');
/*!40000 ALTER TABLE `population` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `project`
--

DROP TABLE IF EXISTS `project`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `project` (
  `row_id` int(10) unsigned NOT NULL auto_increment,
  `project_id` smallint(5) unsigned NOT NULL default '0',
  `ssid` mediumint(8) unsigned default NULL,
  `name` varchar(255) NOT NULL default '',
  `hierarchy_name` varchar(255) NOT NULL default '',
  `study_id` smallint(5) default NULL,
  `note_id` mediumint(8) unsigned default NULL,
  `changed` datetime NOT NULL default '0000-00-00 00:00:00',
  `latest` tinyint(1) default '0',
  PRIMARY KEY  (`row_id`),
  KEY `project_id` (`project_id`),
  KEY `ssid` (`ssid`),
  KEY `latest` (`latest`),
  KEY `hierarchy_name` (`hierarchy_name`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `project`
--

LOCK TABLES `project` WRITE;
/*!40000 ALTER TABLE `project` DISABLE KEYS */;
INSERT INTO `project` VALUES (1,1,NULL,'G0325 [gex] Wellcome Trust Strategic Award application - HIPS','G0325_gex_Wellcome_Trust_Strategic_Award_application_HIPS',NULL,NULL,'2013-08-09 14:07:53',0),(2,1,2625,'G0325 [gex] Wellcome Trust Strategic Award application - HIPS','G0325_gex_Wellcome_Trust_Strategic_Award_application_HIPS',NULL,NULL,'2013-08-09 14:07:53',0),(3,1,2625,'G0325 [gex] Wellcome Trust Strategic Award application - HIPS','G0325_gex_Wellcome_Trust_Strategic_Award_application_HIPS',1,NULL,'2013-08-09 14:07:54',1);
/*!40000 ALTER TABLE `project` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `sample`
--

DROP TABLE IF EXISTS `sample`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `sample` (
  `row_id` int(10) unsigned NOT NULL auto_increment,
  `sample_id` smallint(5) unsigned NOT NULL default '0',
  `project_id` smallint(5) unsigned NOT NULL default '0',
  `ssid` mediumint(8) unsigned default NULL,
  `name` varchar(255) NOT NULL default '',
  `hierarchy_name` varchar(40) NOT NULL default '',
  `individual_id` smallint(5) unsigned default NULL,
  `note_id` mediumint(8) unsigned default NULL,
  `changed` datetime NOT NULL default '0000-00-00 00:00:00',
  `latest` tinyint(1) default '0',
  PRIMARY KEY  (`row_id`),
  KEY `sample_id` (`sample_id`),
  KEY `ssid` (`ssid`),
  KEY `latest` (`latest`),
  KEY `project_id` (`project_id`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `sample`
--

LOCK TABLES `sample` WRITE;
/*!40000 ALTER TABLE `sample` DISABLE KEYS */;
INSERT INTO `sample` VALUES (97,97,0,NULL,'FS11.C','qc1hip5529784',NULL,NULL,'2013-08-09 14:07:57',0),(98,97,1,NULL,'FS11.C','qc1hip5529784',NULL,NULL,'2013-08-09 14:07:57',0),(99,97,1,1625284,'FS11.C','87e7ee6f-e16f-41f6-94c5-194933e2b192',9,2,'2013-08-09 14:07:57',1),(100,100,0,NULL,'FS11.B','qc1hip5529783',NULL,NULL,'2013-08-09 14:07:57',0),(101,100,1,NULL,'FS11.B','qc1hip5529783',NULL,NULL,'2013-08-09 14:07:57',0),(102,100,1,1625283,'FS11.B','3d7a2438-d7bb-40e3-af7b-a0923a253069',9,2,'2013-08-09 14:07:57',1),(103,103,0,NULL,'FS11.A','qc1hip5529782',NULL,NULL,'2013-08-09 14:07:57',0),(104,103,1,NULL,'FS11.A','qc1hip5529782',NULL,NULL,'2013-08-09 14:07:57',0),(105,103,1,1625282,'FS11.A','59e7eb47-c910-4ca8-9cf0-d829539db05b',9,2,'2013-08-09 14:07:57',1),(112,112,0,NULL,'FS11_Control_PBMC','qc1hip5529785',NULL,NULL,'2013-08-09 14:07:58',0),(113,112,1,NULL,'FS11_Control_PBMC','qc1hip5529785',NULL,NULL,'2013-08-09 14:07:58',0),(114,112,1,1625285,'FS11_Control_PBMC','4689a6b5-4dc3-4532-a546-565b41eaa8ba',9,1,'2013-08-09 14:07:58',1);
/*!40000 ALTER TABLE `sample` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `schema_version`
--

DROP TABLE IF EXISTS `schema_version`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `schema_version` (
  `schema_version` mediumint(8) unsigned NOT NULL,
  PRIMARY KEY  (`schema_version`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `schema_version`
--

LOCK TABLES `schema_version` WRITE;
/*!40000 ALTER TABLE `schema_version` DISABLE KEYS */;
INSERT INTO `schema_version` VALUES (24);
/*!40000 ALTER TABLE `schema_version` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `seq_centre`
--

DROP TABLE IF EXISTS `seq_centre`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `seq_centre` (
  `seq_centre_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`seq_centre_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `seq_centre`
--

LOCK TABLES `seq_centre` WRITE;
/*!40000 ALTER TABLE `seq_centre` DISABLE KEYS */;
INSERT INTO `seq_centre` VALUES (1,'SC');
/*!40000 ALTER TABLE `seq_centre` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `seq_request`
--

DROP TABLE IF EXISTS `seq_request`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `seq_request` (
  `row_id` int(10) unsigned NOT NULL auto_increment,
  `seq_request_id` mediumint(8) unsigned NOT NULL default '0',
  `library_id` smallint(5) unsigned default NULL,
  `multiplex_pool_id` smallint(5) unsigned default NULL,
  `ssid` mediumint(8) unsigned default NULL,
  `seq_type` enum('HiSeq Paired end sequencing','Illumina-A HiSeq Paired end sequencing','Illumina-A Paired end sequencing','Illumina-A Pulldown ISC','Illumina-A Pulldown SC','Illumina-A Pulldown WGS','Illumina-A Single ended hi seq sequencing','Illumina-A Single ended sequencing','Illumina-B HiSeq Paired end sequencing','Illumina-B Paired end sequencing','Illumina-B Single ended hi seq sequencing','Illumina-B Single ended sequencing','Illumina-C HiSeq Paired end sequencing','Illumina-C MiSeq sequencing','Illumina-C Paired end sequencing','Illumina-C Single ended hi seq sequencing','Illumina-C Single ended sequencing','MiSeq sequencing','Paired end sequencing','Single ended hi seq sequencing','Single Ended Hi Seq Sequencing Control','Single ended sequencing') default 'Single ended sequencing',
  `seq_status` enum('unknown','pending','started','passed','failed','cancelled','hold') default 'unknown',
  `note_id` mediumint(8) unsigned default NULL,
  `changed` datetime NOT NULL default '0000-00-00 00:00:00',
  `latest` tinyint(1) default '0',
  PRIMARY KEY  (`row_id`),
  KEY `seq_request_id` (`seq_request_id`),
  KEY `ssid` (`ssid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `seq_request`
--

LOCK TABLES `seq_request` WRITE;
/*!40000 ALTER TABLE `seq_request` DISABLE KEYS */;
/*!40000 ALTER TABLE `seq_request` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `seq_tech`
--

DROP TABLE IF EXISTS `seq_tech`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `seq_tech` (
  `seq_tech_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`seq_tech_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `seq_tech`
--

LOCK TABLES `seq_tech` WRITE;
/*!40000 ALTER TABLE `seq_tech` DISABLE KEYS */;
INSERT INTO `seq_tech` VALUES (1,'SLX');
/*!40000 ALTER TABLE `seq_tech` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `species`
--

DROP TABLE IF EXISTS `species`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `species` (
  `species_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL,
  `taxon_id` mediumint(8) unsigned NOT NULL,
  PRIMARY KEY  (`species_id`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `species`
--

LOCK TABLES `species` WRITE;
/*!40000 ALTER TABLE `species` DISABLE KEYS */;
INSERT INTO `species` VALUES (1,'Homo sapiens',9606);
/*!40000 ALTER TABLE `species` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `study`
--

DROP TABLE IF EXISTS `study`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `study` (
  `study_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL default '',
  `acc` varchar(40) default NULL,
  `ssid` mediumint(8) unsigned default NULL,
  `note_id` mediumint(8) unsigned default NULL,
  PRIMARY KEY  (`study_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `study`
--

LOCK TABLES `study` WRITE;
/*!40000 ALTER TABLE `study` DISABLE KEYS */;
INSERT INTO `study` VALUES (1,'','2625',NULL,NULL);
/*!40000 ALTER TABLE `study` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `submission`
--

DROP TABLE IF EXISTS `submission`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `submission` (
  `submission_id` smallint(5) unsigned NOT NULL auto_increment,
  `date` datetime NOT NULL default '0000-00-00 00:00:00',
  `name` varchar(255) NOT NULL default '',
  `acc` varchar(40) default NULL,
  PRIMARY KEY  (`submission_id`),
  UNIQUE KEY `acc` (`acc`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `submission`
--

LOCK TABLES `submission` WRITE;
/*!40000 ALTER TABLE `submission` DISABLE KEYS */;
/*!40000 ALTER TABLE `submission` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Final view structure for view `latest_file`
--

/*!50001 DROP TABLE `latest_file`*/;
/*!50001 DROP VIEW IF EXISTS `latest_file`*/;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`vreseq_rw`@`%` SQL SECURITY DEFINER */
/*!50001 VIEW `latest_file` AS select `file`.`row_id` AS `row_id`,`file`.`file_id` AS `file_id`,`file`.`lane_id` AS `lane_id`,`file`.`name` AS `name`,`file`.`hierarchy_name` AS `hierarchy_name`,`file`.`processed` AS `processed`,`file`.`type` AS `type`,`file`.`readlen` AS `readlen`,`file`.`raw_reads` AS `raw_reads`,`file`.`raw_bases` AS `raw_bases`,`file`.`mean_q` AS `mean_q`,`file`.`md5` AS `md5`,`file`.`note_id` AS `note_id`,`file`.`changed` AS `changed`,`file`.`latest` AS `latest` from `file` where (`file`.`latest` = 1) */;

--
-- Final view structure for view `latest_lane`
--

/*!50001 DROP TABLE `latest_lane`*/;
/*!50001 DROP VIEW IF EXISTS `latest_lane`*/;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`vreseq_rw`@`%` SQL SECURITY DEFINER */
/*!50001 VIEW `latest_lane` AS select `lane`.`row_id` AS `row_id`,`lane`.`lane_id` AS `lane_id`,`lane`.`library_id` AS `library_id`,`lane`.`seq_request_id` AS `seq_request_id`,`lane`.`name` AS `name`,`lane`.`hierarchy_name` AS `hierarchy_name`,`lane`.`acc` AS `acc`,`lane`.`readlen` AS `readlen`,`lane`.`paired` AS `paired`,`lane`.`raw_reads` AS `raw_reads`,`lane`.`raw_bases` AS `raw_bases`,`lane`.`npg_qc_status` AS `npg_qc_status`,`lane`.`processed` AS `processed`,`lane`.`auto_qc_status` AS `auto_qc_status`,`lane`.`qc_status` AS `qc_status`,`lane`.`gt_status` AS `gt_status`,`lane`.`submission_id` AS `submission_id`,`lane`.`withdrawn` AS `withdrawn`,`lane`.`note_id` AS `note_id`,`lane`.`changed` AS `changed`,`lane`.`run_date` AS `run_date`,`lane`.`storage_path` AS `storage_path`,`lane`.`latest` AS `latest` from `lane` where (`lane`.`latest` = 1) */;

--
-- Final view structure for view `latest_library`
--

/*!50001 DROP TABLE `latest_library`*/;
/*!50001 DROP VIEW IF EXISTS `latest_library`*/;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`vreseq_rw`@`%` SQL SECURITY DEFINER */
/*!50001 VIEW `latest_library` AS select `library`.`row_id` AS `row_id`,`library`.`library_id` AS `library_id`,`library`.`library_request_id` AS `library_request_id`,`library`.`sample_id` AS `sample_id`,`library`.`ssid` AS `ssid`,`library`.`name` AS `name`,`library`.`hierarchy_name` AS `hierarchy_name`,`library`.`prep_status` AS `prep_status`,`library`.`auto_qc_status` AS `auto_qc_status`,`library`.`qc_status` AS `qc_status`,`library`.`fragment_size_from` AS `fragment_size_from`,`library`.`fragment_size_to` AS `fragment_size_to`,`library`.`library_type_id` AS `library_type_id`,`library`.`library_tag` AS `library_tag`,`library`.`library_tag_group` AS `library_tag_group`,`library`.`library_tag_sequence` AS `library_tag_sequence`,`library`.`seq_centre_id` AS `seq_centre_id`,`library`.`seq_tech_id` AS `seq_tech_id`,`library`.`open` AS `open`,`library`.`note_id` AS `note_id`,`library`.`changed` AS `changed`,`library`.`latest` AS `latest` from `library` where (`library`.`latest` = 1) */;

--
-- Final view structure for view `latest_library_request`
--

/*!50001 DROP TABLE `latest_library_request`*/;
/*!50001 DROP VIEW IF EXISTS `latest_library_request`*/;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`vreseq_rw`@`%` SQL SECURITY DEFINER */
/*!50001 VIEW `latest_library_request` AS select `library_request`.`row_id` AS `row_id`,`library_request`.`library_request_id` AS `library_request_id`,`library_request`.`sample_id` AS `sample_id`,`library_request`.`ssid` AS `ssid`,`library_request`.`prep_status` AS `prep_status`,`library_request`.`note_id` AS `note_id`,`library_request`.`changed` AS `changed`,`library_request`.`latest` AS `latest` from `library_request` where (`library_request`.`latest` = 1) */;

--
-- Final view structure for view `latest_mapstats`
--

/*!50001 DROP TABLE `latest_mapstats`*/;
/*!50001 DROP VIEW IF EXISTS `latest_mapstats`*/;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`vreseq_rw`@`%` SQL SECURITY DEFINER */
/*!50001 VIEW `latest_mapstats` AS select `mapstats`.`row_id` AS `row_id`,`mapstats`.`mapstats_id` AS `mapstats_id`,`mapstats`.`lane_id` AS `lane_id`,`mapstats`.`mapper_id` AS `mapper_id`,`mapstats`.`assembly_id` AS `assembly_id`,`mapstats`.`raw_reads` AS `raw_reads`,`mapstats`.`raw_bases` AS `raw_bases`,`mapstats`.`clip_bases` AS `clip_bases`,`mapstats`.`reads_mapped` AS `reads_mapped`,`mapstats`.`reads_paired` AS `reads_paired`,`mapstats`.`bases_mapped` AS `bases_mapped`,`mapstats`.`rmdup_reads_mapped` AS `rmdup_reads_mapped`,`mapstats`.`rmdup_bases_mapped` AS `rmdup_bases_mapped`,`mapstats`.`adapter_reads` AS `adapter_reads`,`mapstats`.`error_rate` AS `error_rate`,`mapstats`.`mean_insert` AS `mean_insert`,`mapstats`.`sd_insert` AS `sd_insert`,`mapstats`.`gt_expected` AS `gt_expected`,`mapstats`.`gt_found` AS `gt_found`,`mapstats`.`gt_ratio` AS `gt_ratio`,`mapstats`.`note_id` AS `note_id`,`mapstats`.`changed` AS `changed`,`mapstats`.`latest` AS `latest`,`mapstats`.`bait_near_bases_mapped` AS `bait_near_bases_mapped`,`mapstats`.`target_near_bases_mapped` AS `target_near_bases_mapped`,`mapstats`.`bait_bases_mapped` AS `bait_bases_mapped`,`mapstats`.`mean_bait_coverage` AS `mean_bait_coverage`,`mapstats`.`bait_coverage_sd` AS `bait_coverage_sd`,`mapstats`.`off_bait_bases` AS `off_bait_bases`,`mapstats`.`reads_on_bait` AS `reads_on_bait`,`mapstats`.`reads_on_bait_near` AS `reads_on_bait_near`,`mapstats`.`reads_on_target` AS `reads_on_target`,`mapstats`.`reads_on_target_near` AS `reads_on_target_near`,`mapstats`.`target_bases_mapped` AS `target_bases_mapped`,`mapstats`.`mean_target_coverage` AS `mean_target_coverage`,`mapstats`.`target_coverage_sd` AS `target_coverage_sd`,`mapstats`.`target_bases_1X` AS `target_bases_1X`,`mapstats`.`target_bases_2X` AS `target_bases_2X`,`mapstats`.`target_bases_5X` AS `target_bases_5X`,`mapstats`.`target_bases_10X` AS `target_bases_10X`,`mapstats`.`target_bases_20X` AS `target_bases_20X`,`mapstats`.`target_bases_50X` AS `target_bases_50X`,`mapstats`.`target_bases_100X` AS `target_bases_100X`,`mapstats`.`exome_design_id` AS `exome_design_id`,`mapstats`.`percentage_reads_with_transposon` AS `percentage_reads_with_transposon`,`mapstats`.`is_qc` AS `is_qc`,`mapstats`.`prefix` AS `prefix` from `mapstats` where (`mapstats`.`latest` = 1) */;

--
-- Final view structure for view `latest_project`
--

/*!50001 DROP TABLE `latest_project`*/;
/*!50001 DROP VIEW IF EXISTS `latest_project`*/;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`vreseq_rw`@`%` SQL SECURITY DEFINER */
/*!50001 VIEW `latest_project` AS select `project`.`row_id` AS `row_id`,`project`.`project_id` AS `project_id`,`project`.`ssid` AS `ssid`,`project`.`name` AS `name`,`project`.`hierarchy_name` AS `hierarchy_name`,`project`.`study_id` AS `study_id`,`project`.`note_id` AS `note_id`,`project`.`changed` AS `changed`,`project`.`latest` AS `latest` from `project` where (`project`.`latest` = 1) */;

--
-- Final view structure for view `latest_sample`
--

/*!50001 DROP TABLE `latest_sample`*/;
/*!50001 DROP VIEW IF EXISTS `latest_sample`*/;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`vreseq_rw`@`%` SQL SECURITY DEFINER */
/*!50001 VIEW `latest_sample` AS select `sample`.`row_id` AS `row_id`,`sample`.`sample_id` AS `sample_id`,`sample`.`project_id` AS `project_id`,`sample`.`ssid` AS `ssid`,`sample`.`name` AS `name`,`sample`.`hierarchy_name` AS `hierarchy_name`,`sample`.`individual_id` AS `individual_id`,`sample`.`note_id` AS `note_id`,`sample`.`changed` AS `changed`,`sample`.`latest` AS `latest` from `sample` where (`sample`.`latest` = 1) */;

--
-- Final view structure for view `latest_seq_request`
--

/*!50001 DROP TABLE `latest_seq_request`*/;
/*!50001 DROP VIEW IF EXISTS `latest_seq_request`*/;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`vreseq_rw`@`%` SQL SECURITY DEFINER */
/*!50001 VIEW `latest_seq_request` AS select `seq_request`.`row_id` AS `row_id`,`seq_request`.`seq_request_id` AS `seq_request_id`,`seq_request`.`library_id` AS `library_id`,`seq_request`.`multiplex_pool_id` AS `multiplex_pool_id`,`seq_request`.`ssid` AS `ssid`,`seq_request`.`seq_type` AS `seq_type`,`seq_request`.`seq_status` AS `seq_status`,`seq_request`.`note_id` AS `note_id`,`seq_request`.`changed` AS `changed`,`seq_request`.`latest` AS `latest` from `seq_request` where (`seq_request`.`latest` = 1) */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2013-09-02 16:14:58
