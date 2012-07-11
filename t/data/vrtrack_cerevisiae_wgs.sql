-- MySQL dump 10.11
--
-- Host: mcs4a    Database: vrtrack_testdb_for_vrpipe_sb10
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
INSERT INTO `file` VALUES (1,1,0,'7369_5#31.bam','7369_5_31_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(2,1,1,'7369_5#31.bam','7369_5_31.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(3,1,1,'7369_5#31.bam','7369_5_31.bam',0,4,NULL,NULL,NULL,NULL,'29046bc1c9685bba5e9b96a69282d5e0',NULL,'2012-05-02 12:44:06',1),(4,4,0,'7369_5#30.bam','7369_5_30_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(5,4,4,'7369_5#30.bam','7369_5_30.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(6,4,4,'7369_5#30.bam','7369_5_30.bam',0,4,NULL,NULL,NULL,NULL,'21259f7419fca441d7313c7b372addc0',NULL,'2012-05-02 12:44:06',1),(7,7,0,'7369_5#29.bam','7369_5_29_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(8,7,7,'7369_5#29.bam','7369_5_29.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(9,7,7,'7369_5#29.bam','7369_5_29.bam',0,4,NULL,NULL,NULL,NULL,'2e848f50c7986a4932575beab2fabf4d',NULL,'2012-05-02 12:44:06',1),(10,10,0,'7369_5#28.bam','7369_5_28_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(11,10,10,'7369_5#28.bam','7369_5_28.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(12,10,10,'7369_5#28.bam','7369_5_28.bam',0,4,NULL,NULL,NULL,NULL,'ac9e85c2ac9c2782bb6e19ada462d5e0',NULL,'2012-05-02 12:44:06',1),(13,13,0,'7369_5#27.bam','7369_5_27_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(14,13,13,'7369_5#27.bam','7369_5_27.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(15,13,13,'7369_5#27.bam','7369_5_27.bam',0,4,NULL,NULL,NULL,NULL,'76ee5f08e761fa5baf637a6ad383ad19',NULL,'2012-05-02 12:44:06',1),(16,16,0,'7369_5#26.bam','7369_5_26_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(17,16,16,'7369_5#26.bam','7369_5_26.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:06',0),(18,16,16,'7369_5#26.bam','7369_5_26.bam',0,4,NULL,NULL,NULL,NULL,'00be8b84c1d8ac72a094a7d0705db9eb',NULL,'2012-05-02 12:44:06',1),(19,19,0,'7369_5#25.bam','7369_5_25_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(20,19,19,'7369_5#25.bam','7369_5_25.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(21,19,19,'7369_5#25.bam','7369_5_25.bam',0,4,NULL,NULL,NULL,NULL,'6f4192bca4b15f4b22406e2c629c7223',NULL,'2012-05-02 12:44:07',1),(22,22,0,'7369_5#24.bam','7369_5_24_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(23,22,22,'7369_5#24.bam','7369_5_24.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(24,22,22,'7369_5#24.bam','7369_5_24.bam',0,4,NULL,NULL,NULL,NULL,'e689d397d809d0d7d96e6b34a1f8ff08',NULL,'2012-05-02 12:44:07',1),(25,25,0,'7369_5#23.bam','7369_5_23_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(26,25,25,'7369_5#23.bam','7369_5_23.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(27,25,25,'7369_5#23.bam','7369_5_23.bam',0,4,NULL,NULL,NULL,NULL,'1093582a1df138e4dc598ba7b5b3090d',NULL,'2012-05-02 12:44:07',1),(28,28,0,'7369_5#22.bam','7369_5_22_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(29,28,28,'7369_5#22.bam','7369_5_22.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(30,28,28,'7369_5#22.bam','7369_5_22.bam',0,4,NULL,NULL,NULL,NULL,'1e38fc82f1abda844dba66a759e3b169',NULL,'2012-05-02 12:44:07',1),(31,31,0,'7369_5#21.bam','7369_5_21_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(32,31,31,'7369_5#21.bam','7369_5_21.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(33,31,31,'7369_5#21.bam','7369_5_21.bam',0,4,NULL,NULL,NULL,NULL,'44ca95b23e5222cb01d1db540e14cda8',NULL,'2012-05-02 12:44:07',1),(34,34,0,'7369_5#20.bam','7369_5_20_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(35,34,34,'7369_5#20.bam','7369_5_20.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(36,34,34,'7369_5#20.bam','7369_5_20.bam',0,4,NULL,NULL,NULL,NULL,'062630dfaddb976f18db49d856936d5f',NULL,'2012-05-02 12:44:07',1),(37,37,0,'7369_5#18.bam','7369_5_18_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(38,37,37,'7369_5#18.bam','7369_5_18.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(39,37,37,'7369_5#18.bam','7369_5_18.bam',0,4,NULL,NULL,NULL,NULL,'c301f0ee223abbaa44a92fba5b9a2f18',NULL,'2012-05-02 12:44:07',1),(40,40,0,'7369_5#17.bam','7369_5_17_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(41,40,40,'7369_5#17.bam','7369_5_17.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(42,40,40,'7369_5#17.bam','7369_5_17.bam',0,4,NULL,NULL,NULL,NULL,'401626c972aa3fe358cb5fd5a0e04df2',NULL,'2012-05-02 12:44:07',1),(43,43,0,'7369_5#16.bam','7369_5_16_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(44,43,43,'7369_5#16.bam','7369_5_16.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(45,43,43,'7369_5#16.bam','7369_5_16.bam',0,4,NULL,NULL,NULL,NULL,'c9bbe45ca0e38e085ac510fbb5e6d08c',NULL,'2012-05-02 12:44:07',1),(46,46,0,'7369_5#15.bam','7369_5_15_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(47,46,46,'7369_5#15.bam','7369_5_15.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:07',0),(48,46,46,'7369_5#15.bam','7369_5_15.bam',0,4,NULL,NULL,NULL,NULL,'4927755b9bca21f3d3cdfe0301c9535f',NULL,'2012-05-02 12:44:07',1),(49,49,0,'7369_5#14.bam','7369_5_14_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(50,49,49,'7369_5#14.bam','7369_5_14.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(51,49,49,'7369_5#14.bam','7369_5_14.bam',0,4,NULL,NULL,NULL,NULL,'2a6f90cab3c8d0c7249094c5769abeb8',NULL,'2012-05-02 12:44:08',1),(52,52,0,'7369_5#13.bam','7369_5_13_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(53,52,52,'7369_5#13.bam','7369_5_13.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(54,52,52,'7369_5#13.bam','7369_5_13.bam',0,4,NULL,NULL,NULL,NULL,'e663b81c962f8a76d2aec1a96bf2ab24',NULL,'2012-05-02 12:44:08',1),(55,55,0,'7369_5#12.bam','7369_5_12_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(56,55,55,'7369_5#12.bam','7369_5_12.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(57,55,55,'7369_5#12.bam','7369_5_12.bam',0,4,NULL,NULL,NULL,NULL,'220d921a19b3697471b95b2b9646e21c',NULL,'2012-05-02 12:44:08',1),(58,58,0,'7369_5#11.bam','7369_5_11_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(59,58,58,'7369_5#11.bam','7369_5_11.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(60,58,58,'7369_5#11.bam','7369_5_11.bam',0,4,NULL,NULL,NULL,NULL,'4a8abcc9030aaac542743c2434437578',NULL,'2012-05-02 12:44:08',1),(61,61,0,'7369_5#9.bam','7369_5_9_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(62,61,61,'7369_5#9.bam','7369_5_9.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(63,61,61,'7369_5#9.bam','7369_5_9.bam',0,4,NULL,NULL,NULL,NULL,'2f8e348a2f1b8d12e360f0f3ab213297',NULL,'2012-05-02 12:44:08',1),(64,64,0,'7369_5#8.bam','7369_5_8_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(65,64,64,'7369_5#8.bam','7369_5_8.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(66,64,64,'7369_5#8.bam','7369_5_8.bam',0,4,NULL,NULL,NULL,NULL,'25965d5e328090b4808d53819729a3ab',NULL,'2012-05-02 12:44:08',1),(67,67,0,'7369_5#7.bam','7369_5_7_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(68,67,67,'7369_5#7.bam','7369_5_7.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(69,67,67,'7369_5#7.bam','7369_5_7.bam',0,4,NULL,NULL,NULL,NULL,'0819da001a86761add58ac8ce45c2f35',NULL,'2012-05-02 12:44:08',1),(70,70,0,'7369_5#6.bam','7369_5_6_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(71,70,70,'7369_5#6.bam','7369_5_6.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(72,70,70,'7369_5#6.bam','7369_5_6.bam',0,4,NULL,NULL,NULL,NULL,'c5b8e147e1e1f6f597122a251117850e',NULL,'2012-05-02 12:44:08',1),(73,73,0,'7369_5#5.bam','7369_5_5_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(74,73,73,'7369_5#5.bam','7369_5_5.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(75,73,73,'7369_5#5.bam','7369_5_5.bam',0,4,NULL,NULL,NULL,NULL,'ae8388db009b5072520c4fefdbe31635',NULL,'2012-05-02 12:44:08',1),(76,76,0,'7369_5#4.bam','7369_5_4_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(77,76,76,'7369_5#4.bam','7369_5_4.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(78,76,76,'7369_5#4.bam','7369_5_4.bam',0,4,NULL,NULL,NULL,NULL,'f1ad1e84aa155482db8dfdaced8a4d7a',NULL,'2012-05-02 12:44:08',1),(79,79,0,'7369_5#3.bam','7369_5_3_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(80,79,79,'7369_5#3.bam','7369_5_3.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(81,79,79,'7369_5#3.bam','7369_5_3.bam',0,4,NULL,NULL,NULL,NULL,'dd1bfa50b6f4ab36221006d67ccddc3f',NULL,'2012-05-02 12:44:08',1),(82,82,0,'7369_5#2.bam','7369_5_2_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(83,82,82,'7369_5#2.bam','7369_5_2.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(84,82,82,'7369_5#2.bam','7369_5_2.bam',0,4,NULL,NULL,NULL,NULL,'1ac3ae2539eae650d39d7510c1e5742e',NULL,'2012-05-02 12:44:08',1),(85,85,0,'7369_5#1.bam','7369_5_1_bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(86,85,85,'7369_5#1.bam','7369_5_1.bam',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2012-05-02 12:44:08',0),(87,85,85,'7369_5#1.bam','7369_5_1.bam',0,4,NULL,NULL,NULL,NULL,'8194f6c33299784d78e8d16fb05eb1c6',NULL,'2012-05-02 12:44:08',1);
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
  `name` varchar(255) default NULL,
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
INSERT INTO `individual` VALUES (1,'SC_MFY5249248','SC_MFY5249248','','unknown','ERS074391',1,1),(2,'SC_MFY5249247','SC_MFY5249247','','unknown','ERS074390',1,1),(3,'SC_MFY5249246','SC_MFY5249246','','unknown','ERS074389',1,1),(4,'SC_MFY5249245','SC_MFY5249245','','unknown','ERS074388',1,1),(5,'SC_MFY5249244','SC_MFY5249244','','unknown','ERS074387',1,1),(6,'SC_MFY5249243','SC_MFY5249243','','unknown','ERS074386',1,1),(7,'SC_MFY5249242','SC_MFY5249242','','unknown','ERS074385',1,1),(8,'SC_MFY5249241','SC_MFY5249241','','unknown','ERS074384',1,1),(9,'SC_MFY5249240','SC_MFY5249240','','unknown','ERS074383',1,1),(10,'SC_MFY5249239','SC_MFY5249239','','unknown','ERS074382',1,1),(11,'SC_MFY5249238','SC_MFY5249238','','unknown','ERS074381',1,1),(12,'SC_MFY5249237','SC_MFY5249237','','unknown','ERS074380',1,1),(13,'SC_MFY5249235','SC_MFY5249235','','unknown','ERS074378',1,1),(14,'SC_MFY5249234','SC_MFY5249234','','unknown','ERS074377',1,1),(15,'SC_MFY5249233','SC_MFY5249233','','unknown','ERS074376',1,1),(16,'SC_MFY5249232','SC_MFY5249232','','unknown','ERS074375',1,1),(17,'SC_MFY5249231','SC_MFY5249231','','unknown','ERS074374',1,1),(18,'SC_MFY5249230','SC_MFY5249230','','unknown','ERS074373',1,1),(19,'SC_MFY5249229','SC_MFY5249229','','unknown','ERS074372',1,1),(20,'SC_MFY5249228','SC_MFY5249228','','unknown','ERS074371',1,1),(21,'SC_MFY5249226','SC_MFY5249226','','unknown','ERS074369',1,1),(22,'SC_MFY5249225','SC_MFY5249225','','unknown','ERS074368',1,1),(23,'SC_MFY5249224','SC_MFY5249224','','unknown','ERS074367',1,1),(24,'SC_MFY5249223','SC_MFY5249223','','unknown','ERS074366',1,1),(25,'SC_MFY5249222','SC_MFY5249222','','unknown','ERS074365',1,1),(26,'SC_MFY5249221','SC_MFY5249221','','unknown','ERS074364',1,1),(27,'SC_MFY5249220','SC_MFY5249220','','unknown','ERS074363',1,1),(28,'SC_MFY5249219','SC_MFY5249219','','unknown','ERS074362',1,1),(29,'SC_MFY5249218','SC_MFY5249218','','unknown','ERS074361',1,1);
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
INSERT INTO `lane` VALUES (1,1,0,0,'7369_5#31','7369_5_31',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(2,1,1,0,'7369_5#31','7369_5_31',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(3,1,1,0,'7369_5#31','7369_5#31',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,1),(4,4,0,0,'7369_5#30','7369_5_30',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(5,4,5,0,'7369_5#30','7369_5_30',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(6,4,5,0,'7369_5#30','7369_5#30',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,1),(7,7,0,0,'7369_5#29','7369_5_29',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(8,7,9,0,'7369_5#29','7369_5_29',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(9,7,9,0,'7369_5#29','7369_5#29',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,1),(10,10,0,0,'7369_5#28','7369_5_28',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(11,10,13,0,'7369_5#28','7369_5_28',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(12,10,13,0,'7369_5#28','7369_5#28',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,1),(13,13,0,0,'7369_5#27','7369_5_27',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(14,13,17,0,'7369_5#27','7369_5_27',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(15,13,17,0,'7369_5#27','7369_5#27',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,1),(16,16,0,0,'7369_5#26','7369_5_26',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(17,16,21,0,'7369_5#26','7369_5_26',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:06',NULL,NULL,0),(18,16,21,0,'7369_5#26','7369_5#26',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:06',NULL,NULL,1),(19,19,0,0,'7369_5#25','7369_5_25',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(20,19,25,0,'7369_5#25','7369_5_25',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(21,19,25,0,'7369_5#25','7369_5#25',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(22,22,0,0,'7369_5#24','7369_5_24',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(23,22,29,0,'7369_5#24','7369_5_24',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(24,22,29,0,'7369_5#24','7369_5#24',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(25,25,0,0,'7369_5#23','7369_5_23',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(26,25,33,0,'7369_5#23','7369_5_23',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(27,25,33,0,'7369_5#23','7369_5#23',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(28,28,0,0,'7369_5#22','7369_5_22',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(29,28,37,0,'7369_5#22','7369_5_22',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(30,28,37,0,'7369_5#22','7369_5#22',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(31,31,0,0,'7369_5#21','7369_5_21',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(32,31,41,0,'7369_5#21','7369_5_21',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(33,31,41,0,'7369_5#21','7369_5#21',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(34,34,0,0,'7369_5#20','7369_5_20',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(35,34,45,0,'7369_5#20','7369_5_20',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(36,34,45,0,'7369_5#20','7369_5#20',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(37,37,0,0,'7369_5#18','7369_5_18',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(38,37,49,0,'7369_5#18','7369_5_18',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(39,37,49,0,'7369_5#18','7369_5#18',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(40,40,0,0,'7369_5#17','7369_5_17',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(41,40,53,0,'7369_5#17','7369_5_17',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(42,40,53,0,'7369_5#17','7369_5#17',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(43,43,0,0,'7369_5#16','7369_5_16',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(44,43,57,0,'7369_5#16','7369_5_16',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(45,43,57,0,'7369_5#16','7369_5#16',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(46,46,0,0,'7369_5#15','7369_5_15',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(47,46,61,0,'7369_5#15','7369_5_15',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:07',NULL,NULL,0),(48,46,61,0,'7369_5#15','7369_5#15',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:07',NULL,NULL,1),(49,49,0,0,'7369_5#14','7369_5_14',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(50,49,65,0,'7369_5#14','7369_5_14',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(51,49,65,0,'7369_5#14','7369_5#14',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(52,52,0,0,'7369_5#13','7369_5_13',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(53,52,69,0,'7369_5#13','7369_5_13',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(54,52,69,0,'7369_5#13','7369_5#13',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(55,55,0,0,'7369_5#12','7369_5_12',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(56,55,73,0,'7369_5#12','7369_5_12',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(57,55,73,0,'7369_5#12','7369_5#12',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(58,58,0,0,'7369_5#11','7369_5_11',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(59,58,77,0,'7369_5#11','7369_5_11',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(60,58,77,0,'7369_5#11','7369_5#11',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(61,61,0,0,'7369_5#9','7369_5_9',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(62,61,81,0,'7369_5#9','7369_5_9',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(63,61,81,0,'7369_5#9','7369_5#9',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(64,64,0,0,'7369_5#8','7369_5_8',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(65,64,85,0,'7369_5#8','7369_5_8',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(66,64,85,0,'7369_5#8','7369_5#8',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(67,67,0,0,'7369_5#7','7369_5_7',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(68,67,89,0,'7369_5#7','7369_5_7',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(69,67,89,0,'7369_5#7','7369_5#7',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(70,70,0,0,'7369_5#6','7369_5_6',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(71,70,93,0,'7369_5#6','7369_5_6',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(72,70,93,0,'7369_5#6','7369_5#6',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(73,73,0,0,'7369_5#5','7369_5_5',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(74,73,97,0,'7369_5#5','7369_5_5',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(75,73,97,0,'7369_5#5','7369_5#5',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(76,76,0,0,'7369_5#4','7369_5_4',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(77,76,101,0,'7369_5#4','7369_5_4',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(78,76,101,0,'7369_5#4','7369_5#4',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(79,79,0,0,'7369_5#3','7369_5_3',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(80,79,105,0,'7369_5#3','7369_5_3',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(81,79,105,0,'7369_5#3','7369_5#3',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(82,82,0,0,'7369_5#2','7369_5_2',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(83,82,109,0,'7369_5#2','7369_5_2',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(84,82,109,0,'7369_5#2','7369_5#2',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1),(85,85,0,0,'7369_5#1','7369_5_1',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(86,85,113,0,'7369_5#1','7369_5_1',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2012-05-02 12:44:08',NULL,NULL,0),(87,85,113,0,'7369_5#1','7369_5#1',NULL,NULL,1,0,0,'pass',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2012-05-02 12:44:08',NULL,NULL,1);
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
  `seq_type` enum('Single ended sequencing','Paired end sequencing','HiSeq Paired end sequencing','MiSeq sequencing','Single ended hi seq sequencing'),
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
  `name` varchar(255) default NULL,
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
INSERT INTO `library` VALUES (1,1,0,0,NULL,'4103636','4103636','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(2,1,0,1,NULL,'4103636','4103636','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(3,1,0,1,4103636,'4103636','4103636','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(4,1,0,1,4103636,'4103636','4103636','unknown','no_qc','no_qc',313,313,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:06',1),(5,5,0,0,NULL,'4103648','4103648','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(6,5,0,4,NULL,'4103648','4103648','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(7,5,0,4,4103648,'4103648','4103648','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(8,5,0,4,4103648,'4103648','4103648','unknown','no_qc','no_qc',314,314,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:06',1),(9,9,0,0,NULL,'4103660','4103660','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(10,9,0,7,NULL,'4103660','4103660','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(11,9,0,7,4103660,'4103660','4103660','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(12,9,0,7,4103660,'4103660','4103660','unknown','no_qc','no_qc',329,329,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:06',1),(13,13,0,0,NULL,'4103672','4103672','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(14,13,0,10,NULL,'4103672','4103672','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(15,13,0,10,4103672,'4103672','4103672','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(16,13,0,10,4103672,'4103672','4103672','unknown','no_qc','no_qc',332,332,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:06',1),(17,17,0,0,NULL,'4103684','4103684','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(18,17,0,13,NULL,'4103684','4103684','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(19,17,0,13,4103684,'4103684','4103684','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(20,17,0,13,4103684,'4103684','4103684','unknown','no_qc','no_qc',320,320,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:06',1),(21,21,0,0,NULL,'4103696','4103696','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(22,21,0,16,NULL,'4103696','4103696','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(23,21,0,16,4103696,'4103696','4103696','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:06',0),(24,21,0,16,4103696,'4103696','4103696','unknown','no_qc','no_qc',335,335,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:06',1),(25,25,0,0,NULL,'4103708','4103708','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(26,25,0,19,NULL,'4103708','4103708','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(27,25,0,19,4103708,'4103708','4103708','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(28,25,0,19,4103708,'4103708','4103708','unknown','no_qc','no_qc',320,320,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(29,29,0,0,NULL,'4103625','4103625','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(30,29,0,22,NULL,'4103625','4103625','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(31,29,0,22,4103625,'4103625','4103625','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(32,29,0,22,4103625,'4103625','4103625','unknown','no_qc','no_qc',340,340,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(33,33,0,0,NULL,'4103637','4103637','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(34,33,0,25,NULL,'4103637','4103637','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(35,33,0,25,4103637,'4103637','4103637','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(36,33,0,25,4103637,'4103637','4103637','unknown','no_qc','no_qc',312,312,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(37,37,0,0,NULL,'4103649','4103649','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(38,37,0,28,NULL,'4103649','4103649','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(39,37,0,28,4103649,'4103649','4103649','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(40,37,0,28,4103649,'4103649','4103649','unknown','no_qc','no_qc',325,325,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(41,41,0,0,NULL,'4103661','4103661','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(42,41,0,31,NULL,'4103661','4103661','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(43,41,0,31,4103661,'4103661','4103661','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(44,41,0,31,4103661,'4103661','4103661','unknown','no_qc','no_qc',323,323,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(45,45,0,0,NULL,'4103673','4103673','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(46,45,0,34,NULL,'4103673','4103673','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(47,45,0,34,4103673,'4103673','4103673','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(48,45,0,34,4103673,'4103673','4103673','unknown','no_qc','no_qc',328,328,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(49,49,0,0,NULL,'4103697','4103697','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(50,49,0,37,NULL,'4103697','4103697','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(51,49,0,37,4103697,'4103697','4103697','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(52,49,0,37,4103697,'4103697','4103697','unknown','no_qc','no_qc',316,316,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(53,53,0,0,NULL,'4103709','4103709','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(54,53,0,40,NULL,'4103709','4103709','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(55,53,0,40,4103709,'4103709','4103709','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(56,53,0,40,4103709,'4103709','4103709','unknown','no_qc','no_qc',314,314,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(57,57,0,0,NULL,'4103626','4103626','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(58,57,0,43,NULL,'4103626','4103626','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(59,57,0,43,4103626,'4103626','4103626','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(60,57,0,43,4103626,'4103626','4103626','unknown','no_qc','no_qc',315,315,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(61,61,0,0,NULL,'4103638','4103638','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(62,61,0,46,NULL,'4103638','4103638','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(63,61,0,46,4103638,'4103638','4103638','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:07',0),(64,61,0,46,4103638,'4103638','4103638','unknown','no_qc','no_qc',328,328,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:07',1),(65,65,0,0,NULL,'4103650','4103650','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(66,65,0,49,NULL,'4103650','4103650','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(67,65,0,49,4103650,'4103650','4103650','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(68,65,0,49,4103650,'4103650','4103650','unknown','no_qc','no_qc',330,330,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(69,69,0,0,NULL,'4103662','4103662','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(70,69,0,52,NULL,'4103662','4103662','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(71,69,0,52,4103662,'4103662','4103662','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(72,69,0,52,4103662,'4103662','4103662','unknown','no_qc','no_qc',326,326,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(73,73,0,0,NULL,'4103674','4103674','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(74,73,0,55,NULL,'4103674','4103674','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(75,73,0,55,4103674,'4103674','4103674','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(76,73,0,55,4103674,'4103674','4103674','unknown','no_qc','no_qc',329,329,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(77,77,0,0,NULL,'4103686','4103686','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(78,77,0,58,NULL,'4103686','4103686','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(79,77,0,58,4103686,'4103686','4103686','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(80,77,0,58,4103686,'4103686','4103686','unknown','no_qc','no_qc',323,323,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(81,81,0,0,NULL,'4103710','4103710','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(82,81,0,61,NULL,'4103710','4103710','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(83,81,0,61,4103710,'4103710','4103710','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(84,81,0,61,4103710,'4103710','4103710','unknown','no_qc','no_qc',309,309,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(85,85,0,0,NULL,'4103627','4103627','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(86,85,0,64,NULL,'4103627','4103627','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(87,85,0,64,4103627,'4103627','4103627','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(88,85,0,64,4103627,'4103627','4103627','unknown','no_qc','no_qc',313,313,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(89,89,0,0,NULL,'4103639','4103639','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(90,89,0,67,NULL,'4103639','4103639','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(91,89,0,67,4103639,'4103639','4103639','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(92,89,0,67,4103639,'4103639','4103639','unknown','no_qc','no_qc',309,309,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(93,93,0,0,NULL,'4103651','4103651','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(94,93,0,70,NULL,'4103651','4103651','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(95,93,0,70,4103651,'4103651','4103651','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(96,93,0,70,4103651,'4103651','4103651','unknown','no_qc','no_qc',332,332,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(97,97,0,0,NULL,'4103663','4103663','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(98,97,0,73,NULL,'4103663','4103663','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(99,97,0,73,4103663,'4103663','4103663','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(100,97,0,73,4103663,'4103663','4103663','unknown','no_qc','no_qc',319,319,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(101,101,0,0,NULL,'4103675','4103675','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(102,101,0,76,NULL,'4103675','4103675','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(103,101,0,76,4103675,'4103675','4103675','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(104,101,0,76,4103675,'4103675','4103675','unknown','no_qc','no_qc',317,317,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(105,105,0,0,NULL,'4103687','4103687','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(106,105,0,79,NULL,'4103687','4103687','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(107,105,0,79,4103687,'4103687','4103687','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(108,105,0,79,4103687,'4103687','4103687','unknown','no_qc','no_qc',319,319,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(109,109,0,0,NULL,'4103699','4103699','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(110,109,0,82,NULL,'4103699','4103699','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(111,109,0,82,4103699,'4103699','4103699','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(112,109,0,82,4103699,'4103699','4103699','unknown','no_qc','no_qc',317,317,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1),(113,113,0,0,NULL,'4103711','4103711','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(114,113,0,85,NULL,'4103711','4103711','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(115,113,0,85,4103711,'4103711','4103711','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2012-05-02 12:44:08',0),(116,113,0,85,4103711,'4103711','4103711','unknown','no_qc','no_qc',317,317,NULL,NULL,NULL,NULL,1,1,1,NULL,'2012-05-02 12:44:08',1);
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
  `name` varchar(255) default NULL,
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
  `name` varchar(255) default NULL,
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
  `name` varchar(255) default NULL,
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
INSERT INTO `project` VALUES (1,1,NULL,'SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants','SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants',NULL,NULL,'2012-05-02 12:42:51',0),(2,1,2038,'SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants','SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants',NULL,NULL,'2012-05-02 12:42:51',0),(3,1,2038,'SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants','SEQCAP_WGS_Identification_of_mutational_spectra_in_fission_yeast_DNA_repair_and_chromatin_mutants',1,NULL,'2012-05-02 12:42:51',1);
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
  `name` varchar(255) default NULL,
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
INSERT INTO `sample` VALUES (1,1,0,NULL,'SC_MFY5249248','SC_MFY5249248',NULL,NULL,'2012-05-02 12:44:06',0),(2,1,1,NULL,'SC_MFY5249248','SC_MFY5249248',NULL,NULL,'2012-05-02 12:44:06',0),(3,1,1,1306197,'SC_MFY5249248','SC_MFY5249248',1,NULL,'2012-05-02 12:44:06',1),(4,4,0,NULL,'SC_MFY5249247','SC_MFY5249247',NULL,NULL,'2012-05-02 12:44:06',0),(5,4,1,NULL,'SC_MFY5249247','SC_MFY5249247',NULL,NULL,'2012-05-02 12:44:06',0),(6,4,1,1306196,'SC_MFY5249247','SC_MFY5249247',2,NULL,'2012-05-02 12:44:06',1),(7,7,0,NULL,'SC_MFY5249246','SC_MFY5249246',NULL,NULL,'2012-05-02 12:44:06',0),(8,7,1,NULL,'SC_MFY5249246','SC_MFY5249246',NULL,NULL,'2012-05-02 12:44:06',0),(9,7,1,1306195,'SC_MFY5249246','SC_MFY5249246',3,NULL,'2012-05-02 12:44:06',1),(10,10,0,NULL,'SC_MFY5249245','SC_MFY5249245',NULL,NULL,'2012-05-02 12:44:06',0),(11,10,1,NULL,'SC_MFY5249245','SC_MFY5249245',NULL,NULL,'2012-05-02 12:44:06',0),(12,10,1,1306194,'SC_MFY5249245','SC_MFY5249245',4,NULL,'2012-05-02 12:44:06',1),(13,13,0,NULL,'SC_MFY5249244','SC_MFY5249244',NULL,NULL,'2012-05-02 12:44:06',0),(14,13,1,NULL,'SC_MFY5249244','SC_MFY5249244',NULL,NULL,'2012-05-02 12:44:06',0),(15,13,1,1306193,'SC_MFY5249244','SC_MFY5249244',5,NULL,'2012-05-02 12:44:06',1),(16,16,0,NULL,'SC_MFY5249243','SC_MFY5249243',NULL,NULL,'2012-05-02 12:44:06',0),(17,16,1,NULL,'SC_MFY5249243','SC_MFY5249243',NULL,NULL,'2012-05-02 12:44:06',0),(18,16,1,1306192,'SC_MFY5249243','SC_MFY5249243',6,NULL,'2012-05-02 12:44:06',1),(19,19,0,NULL,'SC_MFY5249242','SC_MFY5249242',NULL,NULL,'2012-05-02 12:44:07',0),(20,19,1,NULL,'SC_MFY5249242','SC_MFY5249242',NULL,NULL,'2012-05-02 12:44:07',0),(21,19,1,1306191,'SC_MFY5249242','SC_MFY5249242',7,NULL,'2012-05-02 12:44:07',1),(22,22,0,NULL,'SC_MFY5249241','SC_MFY5249241',NULL,NULL,'2012-05-02 12:44:07',0),(23,22,1,NULL,'SC_MFY5249241','SC_MFY5249241',NULL,NULL,'2012-05-02 12:44:07',0),(24,22,1,1306190,'SC_MFY5249241','SC_MFY5249241',8,NULL,'2012-05-02 12:44:07',1),(25,25,0,NULL,'SC_MFY5249240','SC_MFY5249240',NULL,NULL,'2012-05-02 12:44:07',0),(26,25,1,NULL,'SC_MFY5249240','SC_MFY5249240',NULL,NULL,'2012-05-02 12:44:07',0),(27,25,1,1306189,'SC_MFY5249240','SC_MFY5249240',9,NULL,'2012-05-02 12:44:07',1),(28,28,0,NULL,'SC_MFY5249239','SC_MFY5249239',NULL,NULL,'2012-05-02 12:44:07',0),(29,28,1,NULL,'SC_MFY5249239','SC_MFY5249239',NULL,NULL,'2012-05-02 12:44:07',0),(30,28,1,1306188,'SC_MFY5249239','SC_MFY5249239',10,NULL,'2012-05-02 12:44:07',1),(31,31,0,NULL,'SC_MFY5249238','SC_MFY5249238',NULL,NULL,'2012-05-02 12:44:07',0),(32,31,1,NULL,'SC_MFY5249238','SC_MFY5249238',NULL,NULL,'2012-05-02 12:44:07',0),(33,31,1,1306187,'SC_MFY5249238','SC_MFY5249238',11,NULL,'2012-05-02 12:44:07',1),(34,34,0,NULL,'SC_MFY5249237','SC_MFY5249237',NULL,NULL,'2012-05-02 12:44:07',0),(35,34,1,NULL,'SC_MFY5249237','SC_MFY5249237',NULL,NULL,'2012-05-02 12:44:07',0),(36,34,1,1306186,'SC_MFY5249237','SC_MFY5249237',12,NULL,'2012-05-02 12:44:07',1),(37,37,0,NULL,'SC_MFY5249235','SC_MFY5249235',NULL,NULL,'2012-05-02 12:44:07',0),(38,37,1,NULL,'SC_MFY5249235','SC_MFY5249235',NULL,NULL,'2012-05-02 12:44:07',0),(39,37,1,1306184,'SC_MFY5249235','SC_MFY5249235',13,NULL,'2012-05-02 12:44:07',1),(40,40,0,NULL,'SC_MFY5249234','SC_MFY5249234',NULL,NULL,'2012-05-02 12:44:07',0),(41,40,1,NULL,'SC_MFY5249234','SC_MFY5249234',NULL,NULL,'2012-05-02 12:44:07',0),(42,40,1,1306183,'SC_MFY5249234','SC_MFY5249234',14,NULL,'2012-05-02 12:44:07',1),(43,43,0,NULL,'SC_MFY5249233','SC_MFY5249233',NULL,NULL,'2012-05-02 12:44:07',0),(44,43,1,NULL,'SC_MFY5249233','SC_MFY5249233',NULL,NULL,'2012-05-02 12:44:07',0),(45,43,1,1306182,'SC_MFY5249233','SC_MFY5249233',15,NULL,'2012-05-02 12:44:07',1),(46,46,0,NULL,'SC_MFY5249232','SC_MFY5249232',NULL,NULL,'2012-05-02 12:44:07',0),(47,46,1,NULL,'SC_MFY5249232','SC_MFY5249232',NULL,NULL,'2012-05-02 12:44:07',0),(48,46,1,1306181,'SC_MFY5249232','SC_MFY5249232',16,NULL,'2012-05-02 12:44:07',1),(49,49,0,NULL,'SC_MFY5249231','SC_MFY5249231',NULL,NULL,'2012-05-02 12:44:08',0),(50,49,1,NULL,'SC_MFY5249231','SC_MFY5249231',NULL,NULL,'2012-05-02 12:44:08',0),(51,49,1,1306180,'SC_MFY5249231','SC_MFY5249231',17,NULL,'2012-05-02 12:44:08',1),(52,52,0,NULL,'SC_MFY5249230','SC_MFY5249230',NULL,NULL,'2012-05-02 12:44:08',0),(53,52,1,NULL,'SC_MFY5249230','SC_MFY5249230',NULL,NULL,'2012-05-02 12:44:08',0),(54,52,1,1306179,'SC_MFY5249230','SC_MFY5249230',18,NULL,'2012-05-02 12:44:08',1),(55,55,0,NULL,'SC_MFY5249229','SC_MFY5249229',NULL,NULL,'2012-05-02 12:44:08',0),(56,55,1,NULL,'SC_MFY5249229','SC_MFY5249229',NULL,NULL,'2012-05-02 12:44:08',0),(57,55,1,1306178,'SC_MFY5249229','SC_MFY5249229',19,NULL,'2012-05-02 12:44:08',1),(58,58,0,NULL,'SC_MFY5249228','SC_MFY5249228',NULL,NULL,'2012-05-02 12:44:08',0),(59,58,1,NULL,'SC_MFY5249228','SC_MFY5249228',NULL,NULL,'2012-05-02 12:44:08',0),(60,58,1,1306177,'SC_MFY5249228','SC_MFY5249228',20,NULL,'2012-05-02 12:44:08',1),(61,61,0,NULL,'SC_MFY5249226','SC_MFY5249226',NULL,NULL,'2012-05-02 12:44:08',0),(62,61,1,NULL,'SC_MFY5249226','SC_MFY5249226',NULL,NULL,'2012-05-02 12:44:08',0),(63,61,1,1306175,'SC_MFY5249226','SC_MFY5249226',21,NULL,'2012-05-02 12:44:08',1),(64,64,0,NULL,'SC_MFY5249225','SC_MFY5249225',NULL,NULL,'2012-05-02 12:44:08',0),(65,64,1,NULL,'SC_MFY5249225','SC_MFY5249225',NULL,NULL,'2012-05-02 12:44:08',0),(66,64,1,1306174,'SC_MFY5249225','SC_MFY5249225',22,NULL,'2012-05-02 12:44:08',1),(67,67,0,NULL,'SC_MFY5249224','SC_MFY5249224',NULL,NULL,'2012-05-02 12:44:08',0),(68,67,1,NULL,'SC_MFY5249224','SC_MFY5249224',NULL,NULL,'2012-05-02 12:44:08',0),(69,67,1,1306173,'SC_MFY5249224','SC_MFY5249224',23,NULL,'2012-05-02 12:44:08',1),(70,70,0,NULL,'SC_MFY5249223','SC_MFY5249223',NULL,NULL,'2012-05-02 12:44:08',0),(71,70,1,NULL,'SC_MFY5249223','SC_MFY5249223',NULL,NULL,'2012-05-02 12:44:08',0),(72,70,1,1306172,'SC_MFY5249223','SC_MFY5249223',24,NULL,'2012-05-02 12:44:08',1),(73,73,0,NULL,'SC_MFY5249222','SC_MFY5249222',NULL,NULL,'2012-05-02 12:44:08',0),(74,73,1,NULL,'SC_MFY5249222','SC_MFY5249222',NULL,NULL,'2012-05-02 12:44:08',0),(75,73,1,1306171,'SC_MFY5249222','SC_MFY5249222',25,NULL,'2012-05-02 12:44:08',1),(76,76,0,NULL,'SC_MFY5249221','SC_MFY5249221',NULL,NULL,'2012-05-02 12:44:08',0),(77,76,1,NULL,'SC_MFY5249221','SC_MFY5249221',NULL,NULL,'2012-05-02 12:44:08',0),(78,76,1,1306170,'SC_MFY5249221','SC_MFY5249221',26,NULL,'2012-05-02 12:44:08',1),(79,79,0,NULL,'SC_MFY5249220','SC_MFY5249220',NULL,NULL,'2012-05-02 12:44:08',0),(80,79,1,NULL,'SC_MFY5249220','SC_MFY5249220',NULL,NULL,'2012-05-02 12:44:08',0),(81,79,1,1306169,'SC_MFY5249220','SC_MFY5249220',27,NULL,'2012-05-02 12:44:08',1),(82,82,0,NULL,'SC_MFY5249219','SC_MFY5249219',NULL,NULL,'2012-05-02 12:44:08',0),(83,82,1,NULL,'SC_MFY5249219','SC_MFY5249219',NULL,NULL,'2012-05-02 12:44:08',0),(84,82,1,1306168,'SC_MFY5249219','SC_MFY5249219',28,NULL,'2012-05-02 12:44:08',1),(85,85,0,NULL,'SC_MFY5249218','SC_MFY5249218',NULL,NULL,'2012-05-02 12:44:08',0),(86,85,1,NULL,'SC_MFY5249218','SC_MFY5249218',NULL,NULL,'2012-05-02 12:44:08',0),(87,85,1,1306167,'SC_MFY5249218','SC_MFY5249218',29,NULL,'2012-05-02 12:44:08',1);
/*!40000 ALTER TABLE `sample` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `schema_version`
--

DROP TABLE IF EXISTS `schema_version`;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
CREATE TABLE `schema_version` (
  `schema_version` mediumint(8) unsigned NOT NULL default '0',
  PRIMARY KEY  (`schema_version`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `schema_version`
--

LOCK TABLES `schema_version` WRITE;
/*!40000 ALTER TABLE `schema_version` DISABLE KEYS */;
INSERT INTO `schema_version` VALUES (20);
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
  `name` varchar(255) default NULL,
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
  `seq_type` enum('Single ended sequencing','Paired end sequencing','HiSeq Paired end sequencing','MiSeq sequencing','Single ended hi seq sequencing') default 'Single ended sequencing',
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
  `name` varchar(255) default NULL,
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
  `name` varchar(255) NOT NULL default '',
  `taxon_id` mediumint(8) unsigned NOT NULL default '0',
  PRIMARY KEY  (`species_id`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
SET character_set_client = @saved_cs_client;

--
-- Dumping data for table `species`
--

LOCK TABLES `species` WRITE;
/*!40000 ALTER TABLE `species` DISABLE KEYS */;
INSERT INTO `species` VALUES (1,'Pombe',1);
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
  `name` varchar(255) default NULL,
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
INSERT INTO `study` VALUES (1,'','ERP001017',NULL,NULL);
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
  `name` varchar(255) default NULL,
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

-- Dump completed on 2012-07-11 14:19:10
