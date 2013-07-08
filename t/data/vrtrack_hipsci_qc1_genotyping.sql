-- MySQL dump 10.11
--
-- Host: mcs10    Database: vrtrack_hipsci_qc1_genotyping
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
INSERT INTO `file` VALUES (1,1,0,'9300870166_R06C01.gtc','9300870166_R06C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:49',0),(2,1,1,'9300870166_R06C01.gtc','9300870166_R06C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:49',0),(3,1,1,'9300870166_R06C01.gtc','9300870166_R06C01.gtc',0,7,NULL,NULL,NULL,NULL,'7793f115dadaa5e0a2b4aae5aca89ce9',NULL,'2013-07-01 11:17:49',1),(4,4,0,'9300870166_R05C01.gtc','9300870166_R05C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:49',0),(5,4,4,'9300870166_R05C01.gtc','9300870166_R05C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:49',0),(6,4,4,'9300870166_R05C01.gtc','9300870166_R05C01.gtc',0,7,NULL,NULL,NULL,NULL,'f067b4a2e68db54814afc03d4f99e0bf',NULL,'2013-07-01 11:17:49',1),(7,7,0,'9300870166_R04C01.gtc','9300870166_R04C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(8,7,7,'9300870166_R04C01.gtc','9300870166_R04C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(9,7,7,'9300870166_R04C01.gtc','9300870166_R04C01.gtc',0,7,NULL,NULL,NULL,NULL,'abf124efcb8990f5f1663e0e9eff0b5f',NULL,'2013-07-01 11:17:50',1),(10,10,0,'9300870166_R03C01.gtc','9300870166_R03C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(11,10,10,'9300870166_R03C01.gtc','9300870166_R03C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(12,10,10,'9300870166_R03C01.gtc','9300870166_R03C01.gtc',0,7,NULL,NULL,NULL,NULL,'356e597598e8186332f55828528edc53',NULL,'2013-07-01 11:17:50',1),(13,13,0,'9300870166_R02C02.gtc','9300870166_R02C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(14,13,13,'9300870166_R02C02.gtc','9300870166_R02C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(15,13,13,'9300870166_R02C02.gtc','9300870166_R02C02.gtc',0,7,NULL,NULL,NULL,NULL,'da8c957df59357e591027544a08478b0',NULL,'2013-07-01 11:17:50',1),(16,16,0,'9300870166_R02C01.gtc','9300870166_R02C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(17,16,16,'9300870166_R02C01.gtc','9300870166_R02C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(18,16,16,'9300870166_R02C01.gtc','9300870166_R02C01.gtc',0,7,NULL,NULL,NULL,NULL,'46a6811e1a88e1705c14843b1c7bf41b',NULL,'2013-07-01 11:17:50',1),(19,19,0,'9300870166_R01C02.gtc','9300870166_R01C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(20,19,19,'9300870166_R01C02.gtc','9300870166_R01C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(21,19,19,'9300870166_R01C02.gtc','9300870166_R01C02.gtc',0,7,NULL,NULL,NULL,NULL,'3f483eb779e63351d9fcef53622d3b26',NULL,'2013-07-01 11:17:50',1),(22,22,0,'9300870166_R01C01.gtc','9300870166_R01C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(23,22,22,'9300870166_R01C01.gtc','9300870166_R01C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(24,22,22,'9300870166_R01C01.gtc','9300870166_R01C01.gtc',0,7,NULL,NULL,NULL,NULL,'5b0a0701696b69f0b421d2c9e41a214e',NULL,'2013-07-01 11:17:50',1),(25,25,0,'9300870057_R06C02.gtc','9300870057_R06C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(26,25,25,'9300870057_R06C02.gtc','9300870057_R06C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:50',0),(27,25,25,'9300870057_R06C02.gtc','9300870057_R06C02.gtc',0,7,NULL,NULL,NULL,NULL,'a965097662c37a2d74a7ca429f58f726',NULL,'2013-07-01 11:17:50',1),(28,28,0,'9300870057_R06C01.gtc','9300870057_R06C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(29,28,28,'9300870057_R06C01.gtc','9300870057_R06C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(30,28,28,'9300870057_R06C01.gtc','9300870057_R06C01.gtc',0,7,NULL,NULL,NULL,NULL,'17b7159554bca4ff4376384b385da51f',NULL,'2013-07-01 11:17:51',1),(31,31,0,'9300870057_R05C02.gtc','9300870057_R05C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(32,31,31,'9300870057_R05C02.gtc','9300870057_R05C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(33,31,31,'9300870057_R05C02.gtc','9300870057_R05C02.gtc',0,7,NULL,NULL,NULL,NULL,'c2e7dd2478e9162e07b68b54f9c3dadf',NULL,'2013-07-01 11:17:51',1),(34,34,0,'9300870057_R05C01.gtc','9300870057_R05C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(35,34,34,'9300870057_R05C01.gtc','9300870057_R05C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(36,34,34,'9300870057_R05C01.gtc','9300870057_R05C01.gtc',0,7,NULL,NULL,NULL,NULL,'72742ad6e230b8f5cab46ac489b257a9',NULL,'2013-07-01 11:17:51',1),(37,37,0,'9300870057_R04C02.gtc','9300870057_R04C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(38,37,37,'9300870057_R04C02.gtc','9300870057_R04C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(39,37,37,'9300870057_R04C02.gtc','9300870057_R04C02.gtc',0,7,NULL,NULL,NULL,NULL,'ce61acbe5bb982cd889174914c03fc84',NULL,'2013-07-01 11:17:51',1),(40,40,0,'9300870057_R04C01.gtc','9300870057_R04C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(41,40,40,'9300870057_R04C01.gtc','9300870057_R04C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(42,40,40,'9300870057_R04C01.gtc','9300870057_R04C01.gtc',0,7,NULL,NULL,NULL,NULL,'d8999f9dd5274ec878f0f4f065aadb6a',NULL,'2013-07-01 11:17:51',1),(43,43,0,'9300870057_R03C02.gtc','9300870057_R03C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(44,43,43,'9300870057_R03C02.gtc','9300870057_R03C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(45,43,43,'9300870057_R03C02.gtc','9300870057_R03C02.gtc',0,7,NULL,NULL,NULL,NULL,'3ae47d2ab0641da7215ba2a1ae544a24',NULL,'2013-07-01 11:17:51',1),(46,46,0,'9300870057_R03C01.gtc','9300870057_R03C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(47,46,46,'9300870057_R03C01.gtc','9300870057_R03C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(48,46,46,'9300870057_R03C01.gtc','9300870057_R03C01.gtc',0,7,NULL,NULL,NULL,NULL,'2742c0dfe226b30de47691029b1de794',NULL,'2013-07-01 11:17:51',1),(49,49,0,'9300870057_R02C02.gtc','9300870057_R02C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(50,49,49,'9300870057_R02C02.gtc','9300870057_R02C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(51,49,49,'9300870057_R02C02.gtc','9300870057_R02C02.gtc',0,7,NULL,NULL,NULL,NULL,'b348b99e7a6850738a5ba029554ca6c7',NULL,'2013-07-01 11:17:51',1),(52,52,0,'9300870057_R02C01.gtc','9300870057_R02C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(53,52,52,'9300870057_R02C01.gtc','9300870057_R02C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:51',0),(54,52,52,'9300870057_R02C01.gtc','9300870057_R02C01.gtc',0,7,NULL,NULL,NULL,NULL,'e6837902bfa6222a4a9562102c73c24b',NULL,'2013-07-01 11:17:51',1),(55,55,0,'9300870057_R01C02.gtc','9300870057_R01C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(56,55,55,'9300870057_R01C02.gtc','9300870057_R01C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(57,55,55,'9300870057_R01C02.gtc','9300870057_R01C02.gtc',0,7,NULL,NULL,NULL,NULL,'d892e1682325c88921afa60acdba1ccf',NULL,'2013-07-01 11:17:52',1),(58,58,0,'9300870057_R01C01.gtc','9300870057_R01C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(59,58,58,'9300870057_R01C01.gtc','9300870057_R01C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(60,58,58,'9300870057_R01C01.gtc','9300870057_R01C01.gtc',0,7,NULL,NULL,NULL,NULL,'d7e10a49be4e8b1e42fe71bc68e93856',NULL,'2013-07-01 11:17:52',1),(61,61,0,'9300870039_R06C01.gtc','9300870039_R06C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(62,61,61,'9300870039_R06C01.gtc','9300870039_R06C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(63,61,61,'9300870039_R06C01.gtc','9300870039_R06C01.gtc',0,7,NULL,NULL,NULL,NULL,'2963ef7aef9ad18331a9a70725c7ccf8',NULL,'2013-07-01 11:17:52',1),(64,64,0,'9300870039_R05C01.gtc','9300870039_R05C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(65,64,64,'9300870039_R05C01.gtc','9300870039_R05C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(66,64,64,'9300870039_R05C01.gtc','9300870039_R05C01.gtc',0,7,NULL,NULL,NULL,NULL,'d3702601bc7ef4d6f7f031690c3422cc',NULL,'2013-07-01 11:17:52',1),(67,67,0,'9300870039_R04C01.gtc','9300870039_R04C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(68,67,67,'9300870039_R04C01.gtc','9300870039_R04C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(69,67,67,'9300870039_R04C01.gtc','9300870039_R04C01.gtc',0,7,NULL,NULL,NULL,NULL,'b4b575b2a39f5ee981595b08cea9a511',NULL,'2013-07-01 11:17:52',1),(70,70,0,'9300870039_R03C01.gtc','9300870039_R03C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(71,70,70,'9300870039_R03C01.gtc','9300870039_R03C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(72,70,70,'9300870039_R03C01.gtc','9300870039_R03C01.gtc',0,7,NULL,NULL,NULL,NULL,'1a68559b94ceb883efec0a10f892db65',NULL,'2013-07-01 11:17:52',1),(73,73,0,'9300870039_R02C02.gtc','9300870039_R02C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(74,73,73,'9300870039_R02C02.gtc','9300870039_R02C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(75,73,73,'9300870039_R02C02.gtc','9300870039_R02C02.gtc',0,7,NULL,NULL,NULL,NULL,'150c06c1bad319d029d4549c9e72812b',NULL,'2013-07-01 11:17:52',1),(76,76,0,'9300870039_R02C01.gtc','9300870039_R02C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(77,76,76,'9300870039_R02C01.gtc','9300870039_R02C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(78,76,76,'9300870039_R02C01.gtc','9300870039_R02C01.gtc',0,7,NULL,NULL,NULL,NULL,'40ced944534b76dacce6d6d67c12dd18',NULL,'2013-07-01 11:17:52',1),(79,79,0,'9300870039_R01C02.gtc','9300870039_R01C02_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(80,79,79,'9300870039_R01C02.gtc','9300870039_R01C02.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(81,79,79,'9300870039_R01C02.gtc','9300870039_R01C02.gtc',0,7,NULL,NULL,NULL,NULL,'f7d4a615e4e1ead8378ed428285b3351',NULL,'2013-07-01 11:17:52',1),(82,82,0,'9300870039_R01C01.gtc','9300870039_R01C01_gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(83,82,82,'9300870039_R01C01.gtc','9300870039_R01C01.gtc',0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'2013-07-01 11:17:52',0),(84,82,82,'9300870039_R01C01.gtc','9300870039_R01C01.gtc',0,7,NULL,NULL,NULL,NULL,'b0cec616e02682e826ef7bca828567e2',NULL,'2013-07-01 11:17:52',1);
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
INSERT INTO `individual` VALUES (1,'2a39941c-12b2-41bf-92f3-70b88b66a3a4','2a39941c_12b2_41bf_92f3_70b88b66a3a4','','unknown',NULL,1,1),(2,'647d3009-5603-4b07-bf02-6161c8662f46','647d3009_5603_4b07_bf02_6161c8662f46','','unknown',NULL,1,1),(3,'ca04b23b-c5b0-4389-95a3-5c7c8e6d51f2','ca04b23b_c5b0_4389_95a3_5c7c8e6d51f2','','unknown',NULL,1,1),(4,'27af9a9b-01b2-4cb6-acef-ea52d83e3d26','27af9a9b_01b2_4cb6_acef_ea52d83e3d26','','unknown',NULL,1,1),(5,'6d3d2acf-29a5-41a2-8992-1414706a527d','6d3d2acf_29a5_41a2_8992_1414706a527d','','unknown',NULL,1,1),(6,'7ddb7ee7-0d09-42f1-9aeb-b97933ed6ec5','7ddb7ee7_0d09_42f1_9aeb_b97933ed6ec5','','unknown',NULL,1,1),(7,'787523cc-ac21-4e05-b667-d8a13e9b7d90','787523cc_ac21_4e05_b667_d8a13e9b7d90','','unknown',NULL,1,1);
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
INSERT INTO `lane` VALUES (1,1,0,0,'9300870166_R06C01','9300870166_R06C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:49',NULL,NULL,0),(2,1,1,0,'9300870166_R06C01','9300870166_R06C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:49',NULL,NULL,0),(3,1,1,0,'9300870166_R06C01','9300870166_R06C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:49',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(4,4,0,0,'9300870166_R05C01','9300870166_R05C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:49',NULL,NULL,0),(5,4,6,0,'9300870166_R05C01','9300870166_R05C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:49',NULL,NULL,0),(6,4,6,0,'9300870166_R05C01','9300870166_R05C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:49',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(7,7,0,0,'9300870166_R04C01','9300870166_R04C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(8,7,11,0,'9300870166_R04C01','9300870166_R04C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(9,7,11,0,'9300870166_R04C01','9300870166_R04C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(10,10,0,0,'9300870166_R03C01','9300870166_R03C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(11,10,16,0,'9300870166_R03C01','9300870166_R03C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(12,10,16,0,'9300870166_R03C01','9300870166_R03C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(13,13,0,0,'9300870166_R02C02','9300870166_R02C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(14,13,21,0,'9300870166_R02C02','9300870166_R02C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(15,13,21,0,'9300870166_R02C02','9300870166_R02C02','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(16,16,0,0,'9300870166_R02C01','9300870166_R02C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(17,16,26,0,'9300870166_R02C01','9300870166_R02C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(18,16,26,0,'9300870166_R02C01','9300870166_R02C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(19,19,0,0,'9300870166_R01C02','9300870166_R01C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(20,19,31,0,'9300870166_R01C02','9300870166_R01C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(21,19,31,0,'9300870166_R01C02','9300870166_R01C02','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(22,22,0,0,'9300870166_R01C01','9300870166_R01C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(23,22,36,0,'9300870166_R01C01','9300870166_R01C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(24,22,36,0,'9300870166_R01C01','9300870166_R01C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(25,25,0,0,'9300870057_R06C02','9300870057_R06C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(26,25,41,0,'9300870057_R06C02','9300870057_R06C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,NULL,0),(27,25,41,0,'9300870057_R06C02','9300870057_R06C02','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:50',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(28,28,0,0,'9300870057_R06C01','9300870057_R06C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(29,28,46,0,'9300870057_R06C01','9300870057_R06C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(30,28,46,0,'9300870057_R06C01','9300870057_R06C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(31,31,0,0,'9300870057_R05C02','9300870057_R05C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(32,31,51,0,'9300870057_R05C02','9300870057_R05C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(33,31,51,0,'9300870057_R05C02','9300870057_R05C02','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(34,34,0,0,'9300870057_R05C01','9300870057_R05C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(35,34,56,0,'9300870057_R05C01','9300870057_R05C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(36,34,56,0,'9300870057_R05C01','9300870057_R05C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(37,37,0,0,'9300870057_R04C02','9300870057_R04C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(38,37,61,0,'9300870057_R04C02','9300870057_R04C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(39,37,61,0,'9300870057_R04C02','9300870057_R04C02','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(40,40,0,0,'9300870057_R04C01','9300870057_R04C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(41,40,66,0,'9300870057_R04C01','9300870057_R04C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(42,40,66,0,'9300870057_R04C01','9300870057_R04C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(43,43,0,0,'9300870057_R03C02','9300870057_R03C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(44,43,71,0,'9300870057_R03C02','9300870057_R03C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(45,43,71,0,'9300870057_R03C02','9300870057_R03C02','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(46,46,0,0,'9300870057_R03C01','9300870057_R03C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(47,46,76,0,'9300870057_R03C01','9300870057_R03C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(48,46,76,0,'9300870057_R03C01','9300870057_R03C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(49,49,0,0,'9300870057_R02C02','9300870057_R02C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(50,49,81,0,'9300870057_R02C02','9300870057_R02C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(51,49,81,0,'9300870057_R02C02','9300870057_R02C02','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(52,52,0,0,'9300870057_R02C01','9300870057_R02C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(53,52,86,0,'9300870057_R02C01','9300870057_R02C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,NULL,0),(54,52,86,0,'9300870057_R02C01','9300870057_R02C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:51',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(55,55,0,0,'9300870057_R01C02','9300870057_R01C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(56,55,91,0,'9300870057_R01C02','9300870057_R01C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(57,55,91,0,'9300870057_R01C02','9300870057_R01C02','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(58,58,0,0,'9300870057_R01C01','9300870057_R01C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(59,58,96,0,'9300870057_R01C01','9300870057_R01C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(60,58,96,0,'9300870057_R01C01','9300870057_R01C01','12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/12d6fd7e-bfb8-4383-aee6-aa62c8f8fdab_coreex_hips_20130531.fcr.txt.gz',1),(61,61,0,0,'9300870039_R06C01','9300870039_R06C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(62,61,101,0,'9300870039_R06C01','9300870039_R06C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(63,61,101,0,'9300870039_R06C01','9300870039_R06C01','45a53a77-50bc-4062-b9cb-8dfe82e589f2',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/45a53a77-50bc-4062-b9cb-8dfe82e589f2_coreex_hips_20130613.fcr.txt.gz',1),(64,64,0,0,'9300870039_R05C01','9300870039_R05C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(65,64,105,0,'9300870039_R05C01','9300870039_R05C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(66,64,105,0,'9300870039_R05C01','9300870039_R05C01','45a53a77-50bc-4062-b9cb-8dfe82e589f2',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/45a53a77-50bc-4062-b9cb-8dfe82e589f2_coreex_hips_20130613.fcr.txt.gz',1),(67,67,0,0,'9300870039_R04C01','9300870039_R04C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(68,67,109,0,'9300870039_R04C01','9300870039_R04C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(69,67,109,0,'9300870039_R04C01','9300870039_R04C01','45a53a77-50bc-4062-b9cb-8dfe82e589f2',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/45a53a77-50bc-4062-b9cb-8dfe82e589f2_coreex_hips_20130613.fcr.txt.gz',1),(70,70,0,0,'9300870039_R03C01','9300870039_R03C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(71,70,113,0,'9300870039_R03C01','9300870039_R03C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(72,70,113,0,'9300870039_R03C01','9300870039_R03C01','45a53a77-50bc-4062-b9cb-8dfe82e589f2',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/45a53a77-50bc-4062-b9cb-8dfe82e589f2_coreex_hips_20130613.fcr.txt.gz',1),(73,73,0,0,'9300870039_R02C02','9300870039_R02C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(74,73,117,0,'9300870039_R02C02','9300870039_R02C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(75,73,117,0,'9300870039_R02C02','9300870039_R02C02','45a53a77-50bc-4062-b9cb-8dfe82e589f2',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/45a53a77-50bc-4062-b9cb-8dfe82e589f2_coreex_hips_20130613.fcr.txt.gz',1),(76,76,0,0,'9300870039_R02C01','9300870039_R02C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(77,76,121,0,'9300870039_R02C01','9300870039_R02C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(78,76,121,0,'9300870039_R02C01','9300870039_R02C01','45a53a77-50bc-4062-b9cb-8dfe82e589f2',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/45a53a77-50bc-4062-b9cb-8dfe82e589f2_coreex_hips_20130613.fcr.txt.gz',1),(79,79,0,0,'9300870039_R01C02','9300870039_R01C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(80,79,125,0,'9300870039_R01C02','9300870039_R01C02',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(81,79,125,0,'9300870039_R01C02','9300870039_R01C02','45a53a77-50bc-4062-b9cb-8dfe82e589f2',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/45a53a77-50bc-4062-b9cb-8dfe82e589f2_coreex_hips_20130613.fcr.txt.gz',1),(82,82,0,0,'9300870039_R01C01','9300870039_R01C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(83,82,129,0,'9300870039_R01C01','9300870039_R01C01',NULL,NULL,NULL,NULL,NULL,'pending',0,'no_qc','no_qc','unchecked',NULL,NULL,NULL,'2013-07-01 11:17:52',NULL,NULL,0),(84,82,129,0,'9300870039_R01C01','9300870039_R01C01','45a53a77-50bc-4062-b9cb-8dfe82e589f2',NULL,0,0,0,'pending',0,'no_qc','no_qc','unchecked',NULL,1,NULL,'2013-07-01 11:17:52',NULL,'/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/45a53a77-50bc-4062-b9cb-8dfe82e589f2_coreex_hips_20130613.fcr.txt.gz',1);
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
INSERT INTO `library` VALUES (1,1,0,0,NULL,'283163_B03_qc1hip5533830','283163_B03_qc1hip5533830','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:49',0),(2,1,0,1,NULL,'283163_B03_qc1hip5533830','283163_B03_qc1hip5533830','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:49',0),(3,1,0,1,166475,'283163_B03_qc1hip5533830','283163_B03_qc1hip5533830','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:49',0),(4,1,0,1,166475,'283163_B03_qc1hip5533830','283163_B03_qc1hip5533830','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R06C01',1,1,1,NULL,'2013-07-01 11:17:49',0),(5,1,0,1,166475,'283163_B03_qc1hip5533830','283163_B03_qc1hip5533830','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R06C01',1,1,1,NULL,'2013-07-01 11:17:49',1),(6,6,0,0,NULL,'283163_A03_qc1hip5533829','283163_A03_qc1hip5533829','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:49',0),(7,6,0,4,NULL,'283163_A03_qc1hip5533829','283163_A03_qc1hip5533829','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:49',0),(8,6,0,4,166474,'283163_A03_qc1hip5533829','283163_A03_qc1hip5533829','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:49',0),(9,6,0,4,166474,'283163_A03_qc1hip5533829','283163_A03_qc1hip5533829','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R05C01',1,1,1,NULL,'2013-07-01 11:17:49',0),(10,6,0,4,166474,'283163_A03_qc1hip5533829','283163_A03_qc1hip5533829','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R05C01',1,1,1,NULL,'2013-07-01 11:17:49',1),(11,11,0,0,NULL,'283163_H02_qc1hip5533828','283163_H02_qc1hip5533828','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(12,11,0,7,NULL,'283163_H02_qc1hip5533828','283163_H02_qc1hip5533828','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(13,11,0,7,166473,'283163_H02_qc1hip5533828','283163_H02_qc1hip5533828','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(14,11,0,7,166473,'283163_H02_qc1hip5533828','283163_H02_qc1hip5533828','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R04C01',1,1,1,NULL,'2013-07-01 11:17:50',0),(15,11,0,7,166473,'283163_H02_qc1hip5533828','283163_H02_qc1hip5533828','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R04C01',1,1,1,NULL,'2013-07-01 11:17:50',1),(16,16,0,0,NULL,'283163_G02_qc1hip5533827','283163_G02_qc1hip5533827','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(17,16,0,10,NULL,'283163_G02_qc1hip5533827','283163_G02_qc1hip5533827','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(18,16,0,10,166472,'283163_G02_qc1hip5533827','283163_G02_qc1hip5533827','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(19,16,0,10,166472,'283163_G02_qc1hip5533827','283163_G02_qc1hip5533827','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R03C01',1,1,1,NULL,'2013-07-01 11:17:50',0),(20,16,0,10,166472,'283163_G02_qc1hip5533827','283163_G02_qc1hip5533827','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R03C01',1,1,1,NULL,'2013-07-01 11:17:50',1),(21,21,0,0,NULL,'283163_D03_qc1hip5533832','283163_D03_qc1hip5533832','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(22,21,0,13,NULL,'283163_D03_qc1hip5533832','283163_D03_qc1hip5533832','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(23,21,0,13,166477,'283163_D03_qc1hip5533832','283163_D03_qc1hip5533832','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(24,21,0,13,166477,'283163_D03_qc1hip5533832','283163_D03_qc1hip5533832','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C02',1,1,1,NULL,'2013-07-01 11:17:50',0),(25,21,0,13,166477,'283163_D03_qc1hip5533832','283163_D03_qc1hip5533832','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C02',1,1,1,NULL,'2013-07-01 11:17:50',1),(26,26,0,0,NULL,'283163_F02_qc1hip5533826','283163_F02_qc1hip5533826','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(27,26,0,16,NULL,'283163_F02_qc1hip5533826','283163_F02_qc1hip5533826','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(28,26,0,16,166471,'283163_F02_qc1hip5533826','283163_F02_qc1hip5533826','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(29,26,0,16,166471,'283163_F02_qc1hip5533826','283163_F02_qc1hip5533826','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C01',1,1,1,NULL,'2013-07-01 11:17:50',0),(30,26,0,16,166471,'283163_F02_qc1hip5533826','283163_F02_qc1hip5533826','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C01',1,1,1,NULL,'2013-07-01 11:17:50',1),(31,31,0,0,NULL,'283163_C03_qc1hip5533831','283163_C03_qc1hip5533831','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(32,31,0,19,NULL,'283163_C03_qc1hip5533831','283163_C03_qc1hip5533831','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(33,31,0,19,166476,'283163_C03_qc1hip5533831','283163_C03_qc1hip5533831','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(34,31,0,19,166476,'283163_C03_qc1hip5533831','283163_C03_qc1hip5533831','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C02',1,1,1,NULL,'2013-07-01 11:17:50',0),(35,31,0,19,166476,'283163_C03_qc1hip5533831','283163_C03_qc1hip5533831','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C02',1,1,1,NULL,'2013-07-01 11:17:50',1),(36,36,0,0,NULL,'283163_E02_qc1hip5533825','283163_E02_qc1hip5533825','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(37,36,0,22,NULL,'283163_E02_qc1hip5533825','283163_E02_qc1hip5533825','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(38,36,0,22,166470,'283163_E02_qc1hip5533825','283163_E02_qc1hip5533825','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(39,36,0,22,166470,'283163_E02_qc1hip5533825','283163_E02_qc1hip5533825','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C01',1,1,1,NULL,'2013-07-01 11:17:50',0),(40,36,0,22,166470,'283163_E02_qc1hip5533825','283163_E02_qc1hip5533825','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C01',1,1,1,NULL,'2013-07-01 11:17:50',1),(41,41,0,0,NULL,'283163_D02_qc1hip5533824','283163_D02_qc1hip5533824','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(42,41,0,25,NULL,'283163_D02_qc1hip5533824','283163_D02_qc1hip5533824','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(43,41,0,25,57469,'283163_D02_qc1hip5533824','283163_D02_qc1hip5533824','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:50',0),(44,41,0,25,57469,'283163_D02_qc1hip5533824','283163_D02_qc1hip5533824','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R06C02',1,1,1,NULL,'2013-07-01 11:17:50',0),(45,41,0,25,57469,'283163_D02_qc1hip5533824','283163_D02_qc1hip5533824','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R06C02',1,1,1,NULL,'2013-07-01 11:17:50',1),(46,46,0,0,NULL,'283163_F01_qc1hip5529688','283163_F01_qc1hip5529688','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(47,46,0,28,NULL,'283163_F01_qc1hip5529688','283163_F01_qc1hip5529688','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(48,46,0,28,57188,'283163_F01_qc1hip5529688','283163_F01_qc1hip5529688','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(49,46,0,28,57188,'283163_F01_qc1hip5529688','283163_F01_qc1hip5529688','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R06C01',1,1,1,NULL,'2013-07-01 11:17:51',0),(50,46,0,28,57188,'283163_F01_qc1hip5529688','283163_F01_qc1hip5529688','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R06C01',1,1,1,NULL,'2013-07-01 11:17:51',1),(51,51,0,0,NULL,'283163_C02_qc1hip5533823','283163_C02_qc1hip5533823','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(52,51,0,31,NULL,'283163_C02_qc1hip5533823','283163_C02_qc1hip5533823','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(53,51,0,31,57468,'283163_C02_qc1hip5533823','283163_C02_qc1hip5533823','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(54,51,0,31,57468,'283163_C02_qc1hip5533823','283163_C02_qc1hip5533823','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R05C02',1,1,1,NULL,'2013-07-01 11:17:51',0),(55,51,0,31,57468,'283163_C02_qc1hip5533823','283163_C02_qc1hip5533823','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R05C02',1,1,1,NULL,'2013-07-01 11:17:51',1),(56,56,0,0,NULL,'283163_E01_qc1hip5529687','283163_E01_qc1hip5529687','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(57,56,0,34,NULL,'283163_E01_qc1hip5529687','283163_E01_qc1hip5529687','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(58,56,0,34,57187,'283163_E01_qc1hip5529687','283163_E01_qc1hip5529687','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(59,56,0,34,57187,'283163_E01_qc1hip5529687','283163_E01_qc1hip5529687','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R05C01',1,1,1,NULL,'2013-07-01 11:17:51',0),(60,56,0,34,57187,'283163_E01_qc1hip5529687','283163_E01_qc1hip5529687','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R05C01',1,1,1,NULL,'2013-07-01 11:17:51',1),(61,61,0,0,NULL,'283163_B02_qc1hip5533822','283163_B02_qc1hip5533822','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(62,61,0,37,NULL,'283163_B02_qc1hip5533822','283163_B02_qc1hip5533822','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(63,61,0,37,57467,'283163_B02_qc1hip5533822','283163_B02_qc1hip5533822','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(64,61,0,37,57467,'283163_B02_qc1hip5533822','283163_B02_qc1hip5533822','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R04C02',1,1,1,NULL,'2013-07-01 11:17:51',0),(65,61,0,37,57467,'283163_B02_qc1hip5533822','283163_B02_qc1hip5533822','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R04C02',1,1,1,NULL,'2013-07-01 11:17:51',1),(66,66,0,0,NULL,'283163_D01_qc1hip5529686','283163_D01_qc1hip5529686','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(67,66,0,40,NULL,'283163_D01_qc1hip5529686','283163_D01_qc1hip5529686','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(68,66,0,40,57186,'283163_D01_qc1hip5529686','283163_D01_qc1hip5529686','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(69,66,0,40,57186,'283163_D01_qc1hip5529686','283163_D01_qc1hip5529686','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R04C01',1,1,1,NULL,'2013-07-01 11:17:51',0),(70,66,0,40,57186,'283163_D01_qc1hip5529686','283163_D01_qc1hip5529686','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R04C01',1,1,1,NULL,'2013-07-01 11:17:51',1),(71,71,0,0,NULL,'283163_A02_qc1hip5533821','283163_A02_qc1hip5533821','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(72,71,0,43,NULL,'283163_A02_qc1hip5533821','283163_A02_qc1hip5533821','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(73,71,0,43,57466,'283163_A02_qc1hip5533821','283163_A02_qc1hip5533821','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(74,71,0,43,57466,'283163_A02_qc1hip5533821','283163_A02_qc1hip5533821','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R03C02',1,1,1,NULL,'2013-07-01 11:17:51',0),(75,71,0,43,57466,'283163_A02_qc1hip5533821','283163_A02_qc1hip5533821','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R03C02',1,1,1,NULL,'2013-07-01 11:17:51',1),(76,76,0,0,NULL,'283163_C01_qc1hip5529685','283163_C01_qc1hip5529685','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(77,76,0,46,NULL,'283163_C01_qc1hip5529685','283163_C01_qc1hip5529685','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(78,76,0,46,57185,'283163_C01_qc1hip5529685','283163_C01_qc1hip5529685','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(79,76,0,46,57185,'283163_C01_qc1hip5529685','283163_C01_qc1hip5529685','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R03C01',1,1,1,NULL,'2013-07-01 11:17:51',0),(80,76,0,46,57185,'283163_C01_qc1hip5529685','283163_C01_qc1hip5529685','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R03C01',1,1,1,NULL,'2013-07-01 11:17:51',1),(81,81,0,0,NULL,'283163_H01_qc1hip5533833','283163_H01_qc1hip5533833','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(82,81,0,49,NULL,'283163_H01_qc1hip5533833','283163_H01_qc1hip5533833','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(83,81,0,49,57478,'283163_H01_qc1hip5533833','283163_H01_qc1hip5533833','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(84,81,0,49,57478,'283163_H01_qc1hip5533833','283163_H01_qc1hip5533833','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C02',1,1,1,NULL,'2013-07-01 11:17:51',0),(85,81,0,49,57478,'283163_H01_qc1hip5533833','283163_H01_qc1hip5533833','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C02',1,1,1,NULL,'2013-07-01 11:17:51',1),(86,86,0,0,NULL,'283163_B01_qc1hip5529684','283163_B01_qc1hip5529684','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(87,86,0,52,NULL,'283163_B01_qc1hip5529684','283163_B01_qc1hip5529684','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(88,86,0,52,57184,'283163_B01_qc1hip5529684','283163_B01_qc1hip5529684','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(89,86,0,52,57184,'283163_B01_qc1hip5529684','283163_B01_qc1hip5529684','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C01',1,1,1,NULL,'2013-07-01 11:17:51',0),(90,86,0,52,57184,'283163_B01_qc1hip5529684','283163_B01_qc1hip5529684','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C01',1,1,1,NULL,'2013-07-01 11:17:51',1),(91,91,0,0,NULL,'283163_G01_qc1hip5529689','283163_G01_qc1hip5529689','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(92,91,0,55,NULL,'283163_G01_qc1hip5529689','283163_G01_qc1hip5529689','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(93,91,0,55,57189,'283163_G01_qc1hip5529689','283163_G01_qc1hip5529689','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:51',0),(94,91,0,55,57189,'283163_G01_qc1hip5529689','283163_G01_qc1hip5529689','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C02',1,1,1,NULL,'2013-07-01 11:17:52',0),(95,91,0,55,57189,'283163_G01_qc1hip5529689','283163_G01_qc1hip5529689','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C02',1,1,1,NULL,'2013-07-01 11:17:52',1),(96,96,0,0,NULL,'283163_A01_qc1hip5529683','283163_A01_qc1hip5529683','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(97,96,0,58,NULL,'283163_A01_qc1hip5529683','283163_A01_qc1hip5529683','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(98,96,0,58,57183,'283163_A01_qc1hip5529683','283163_A01_qc1hip5529683','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(99,96,0,58,57183,'283163_A01_qc1hip5529683','283163_A01_qc1hip5529683','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C01',1,1,1,NULL,'2013-07-01 11:17:52',0),(100,96,0,58,57183,'283163_A01_qc1hip5529683','283163_A01_qc1hip5529683','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C01',1,1,1,NULL,'2013-07-01 11:17:52',1),(101,101,0,0,NULL,'285708_F01_qc1hip5539494','285708_F01_qc1hip5539494','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(102,101,0,61,NULL,'285708_F01_qc1hip5539494','285708_F01_qc1hip5539494','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(103,101,0,61,39284,'285708_F01_qc1hip5539494','285708_F01_qc1hip5539494','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(104,101,0,61,39284,'285708_F01_qc1hip5539494','285708_F01_qc1hip5539494','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R06C01',1,1,1,NULL,'2013-07-01 11:17:52',1),(105,105,0,0,NULL,'285708_E01_qc1hip5539493','285708_E01_qc1hip5539493','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(106,105,0,64,NULL,'285708_E01_qc1hip5539493','285708_E01_qc1hip5539493','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(107,105,0,64,39283,'285708_E01_qc1hip5539493','285708_E01_qc1hip5539493','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(108,105,0,64,39283,'285708_E01_qc1hip5539493','285708_E01_qc1hip5539493','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R05C01',1,1,1,NULL,'2013-07-01 11:17:52',1),(109,109,0,0,NULL,'285708_D01_qc1hip5539488','285708_D01_qc1hip5539488','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(110,109,0,67,NULL,'285708_D01_qc1hip5539488','285708_D01_qc1hip5539488','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(111,109,0,67,39278,'285708_D01_qc1hip5539488','285708_D01_qc1hip5539488','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(112,109,0,67,39278,'285708_D01_qc1hip5539488','285708_D01_qc1hip5539488','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R04C01',1,1,1,NULL,'2013-07-01 11:17:52',1),(113,113,0,0,NULL,'285708_C01_qc1hip5539487','285708_C01_qc1hip5539487','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(114,113,0,70,NULL,'285708_C01_qc1hip5539487','285708_C01_qc1hip5539487','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(115,113,0,70,39277,'285708_C01_qc1hip5539487','285708_C01_qc1hip5539487','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(116,113,0,70,39277,'285708_C01_qc1hip5539487','285708_C01_qc1hip5539487','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R03C01',1,1,1,NULL,'2013-07-01 11:17:52',1),(117,117,0,0,NULL,'285708_H01_qc1hip5543904','285708_H01_qc1hip5543904','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(118,117,0,73,NULL,'285708_H01_qc1hip5543904','285708_H01_qc1hip5543904','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(119,117,0,73,39801,'285708_H01_qc1hip5543904','285708_H01_qc1hip5543904','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(120,117,0,73,39801,'285708_H01_qc1hip5543904','285708_H01_qc1hip5543904','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C02',1,1,1,NULL,'2013-07-01 11:17:52',1),(121,121,0,0,NULL,'285708_B01_qc1hip5539486','285708_B01_qc1hip5539486','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(122,121,0,76,NULL,'285708_B01_qc1hip5539486','285708_B01_qc1hip5539486','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(123,121,0,76,39276,'285708_B01_qc1hip5539486','285708_B01_qc1hip5539486','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(124,121,0,76,39276,'285708_B01_qc1hip5539486','285708_B01_qc1hip5539486','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R02C01',1,1,1,NULL,'2013-07-01 11:17:52',1),(125,125,0,0,NULL,'285708_G01_qc1hip5543903','285708_G01_qc1hip5543903','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(126,125,0,79,NULL,'285708_G01_qc1hip5543903','285708_G01_qc1hip5543903','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(127,125,0,79,39800,'285708_G01_qc1hip5543903','285708_G01_qc1hip5543903','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(128,125,0,79,39800,'285708_G01_qc1hip5543903','285708_G01_qc1hip5543903','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C02',1,1,1,NULL,'2013-07-01 11:17:52',1),(129,129,0,0,NULL,'285708_A01_qc1hip5539485','285708_A01_qc1hip5539485','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(130,129,0,82,NULL,'285708_A01_qc1hip5539485','285708_A01_qc1hip5539485','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(131,129,0,82,39275,'285708_A01_qc1hip5539485','285708_A01_qc1hip5539485','unknown','no_qc','no_qc',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,1,NULL,'2013-07-01 11:17:52',0),(132,129,0,82,39275,'285708_A01_qc1hip5539485','285708_A01_qc1hip5539485','unknown','no_qc','no_qc',0,0,NULL,NULL,NULL,'R01C01',1,1,1,NULL,'2013-07-01 11:17:52',1);
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
INSERT INTO `project` VALUES (1,1,NULL,'Wellcome Trust Strategic Award application â€“ HIPS','Wellcome_Trust_Strategic_Award_application_HIPS',NULL,NULL,'2013-07-01 11:17:48',0),(2,1,2624,'Wellcome Trust Strategic Award application â€“ HIPS','Wellcome_Trust_Strategic_Award_application_HIPS',NULL,NULL,'2013-07-01 11:17:49',0),(3,1,2624,'Wellcome Trust Strategic Award application â€“ HIPS','Wellcome_Trust_Strategic_Award_application_HIPS',1,NULL,'2013-07-01 11:17:49',1),(4,4,NULL,'Wellcome Trust Strategic Award application – HIPS','Wellcome_Trust_Strategic_Award_application_HIPS',NULL,NULL,'2013-07-01 11:17:52',0),(5,4,2624,'Wellcome Trust Strategic Award application – HIPS','Wellcome_Trust_Strategic_Award_application_HIPS',NULL,NULL,'2013-07-01 11:17:52',0),(6,4,2624,'Wellcome Trust Strategic Award application – HIPS','Wellcome_Trust_Strategic_Award_application_HIPS',1,NULL,'2013-07-01 11:17:52',1);
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
INSERT INTO `sample` VALUES (1,1,0,NULL,'qc1hip5533830','qc1hip5533830',NULL,NULL,'2013-07-01 11:17:49',0),(2,1,1,NULL,'qc1hip5533830','qc1hip5533830',NULL,NULL,'2013-07-01 11:17:49',0),(3,1,1,1629475,'qc1hip5533830','3972c90e-c37c-4a80-a21c-79ea309fce9d',1,2,'2013-07-01 11:17:49',1),(4,4,0,NULL,'qc1hip5533829','qc1hip5533829',NULL,NULL,'2013-07-01 11:17:49',0),(5,4,1,NULL,'qc1hip5533829','qc1hip5533829',NULL,NULL,'2013-07-01 11:17:49',0),(6,4,1,1629474,'qc1hip5533829','b87bc4ee-7841-4601-9287-25983028f83c',1,2,'2013-07-01 11:17:49',1),(7,7,0,NULL,'qc1hip5533828','qc1hip5533828',NULL,NULL,'2013-07-01 11:17:49',0),(8,7,1,NULL,'qc1hip5533828','qc1hip5533828',NULL,NULL,'2013-07-01 11:17:49',0),(9,7,1,1629473,'qc1hip5533828','b43bb5e0-b938-4c86-a1d8-0162cc8e1eca',2,1,'2013-07-01 11:17:50',1),(10,10,0,NULL,'qc1hip5533827','qc1hip5533827',NULL,NULL,'2013-07-01 11:17:50',0),(11,10,1,NULL,'qc1hip5533827','qc1hip5533827',NULL,NULL,'2013-07-01 11:17:50',0),(12,10,1,1629472,'qc1hip5533827','0bc05be6-8a18-4ea4-af6a-8bb5f31f41c8',2,2,'2013-07-01 11:17:50',1),(13,13,0,NULL,'qc1hip5533832','qc1hip5533832',NULL,NULL,'2013-07-01 11:17:50',0),(14,13,1,NULL,'qc1hip5533832','qc1hip5533832',NULL,NULL,'2013-07-01 11:17:50',0),(15,13,1,1629477,'qc1hip5533832','77e3f1fc-63cd-4ee9-9e2c-0a7ff3474bf1',1,1,'2013-07-01 11:17:50',1),(16,16,0,NULL,'qc1hip5533826','qc1hip5533826',NULL,NULL,'2013-07-01 11:17:50',0),(17,16,1,NULL,'qc1hip5533826','qc1hip5533826',NULL,NULL,'2013-07-01 11:17:50',0),(18,16,1,1629471,'qc1hip5533826','b8227d7e-e993-4210-a9c6-c509a9c11ed9',2,2,'2013-07-01 11:17:50',1),(19,19,0,NULL,'qc1hip5533831','qc1hip5533831',NULL,NULL,'2013-07-01 11:17:50',0),(20,19,1,NULL,'qc1hip5533831','qc1hip5533831',NULL,NULL,'2013-07-01 11:17:50',0),(21,19,1,1629476,'qc1hip5533831','ec257e0f-bc58-4761-9fab-41dc00d182de',1,2,'2013-07-01 11:17:50',1),(22,22,0,NULL,'qc1hip5533825','qc1hip5533825',NULL,NULL,'2013-07-01 11:17:50',0),(23,22,1,NULL,'qc1hip5533825','qc1hip5533825',NULL,NULL,'2013-07-01 11:17:50',0),(24,22,1,1629470,'qc1hip5533825','f629dae2-5cba-44c9-a8f8-ee2bb128506c',2,2,'2013-07-01 11:17:50',1),(25,25,0,NULL,'qc1hip5533824','qc1hip5533824',NULL,NULL,'2013-07-01 11:17:50',0),(26,25,1,NULL,'qc1hip5533824','qc1hip5533824',NULL,NULL,'2013-07-01 11:17:50',0),(27,25,1,1629469,'qc1hip5533824','d95b0f29-07b2-4ca6-a2c3-ea16db230537',3,1,'2013-07-01 11:17:50',1),(28,28,0,NULL,'qc1hip5529688','qc1hip5529688',NULL,NULL,'2013-07-01 11:17:50',0),(29,28,1,NULL,'qc1hip5529688','qc1hip5529688',NULL,NULL,'2013-07-01 11:17:51',0),(30,28,1,1625188,'qc1hip5529688','87e7ee6f-e16f-41f6-94c5-194933e2b192',4,2,'2013-07-01 11:17:51',1),(31,31,0,NULL,'qc1hip5533823','qc1hip5533823',NULL,NULL,'2013-07-01 11:17:51',0),(32,31,1,NULL,'qc1hip5533823','qc1hip5533823',NULL,NULL,'2013-07-01 11:17:51',0),(33,31,1,1629468,'qc1hip5533823','73ac7fdb-773d-4200-97ab-712966d2155a',3,2,'2013-07-01 11:17:51',1),(34,34,0,NULL,'qc1hip5529687','qc1hip5529687',NULL,NULL,'2013-07-01 11:17:51',0),(35,34,1,NULL,'qc1hip5529687','qc1hip5529687',NULL,NULL,'2013-07-01 11:17:51',0),(36,34,1,1625187,'qc1hip5529687','3d7a2438-d7bb-40e3-af7b-a0923a253069',4,2,'2013-07-01 11:17:51',1),(37,37,0,NULL,'qc1hip5533822','qc1hip5533822',NULL,NULL,'2013-07-01 11:17:51',0),(38,37,1,NULL,'qc1hip5533822','qc1hip5533822',NULL,NULL,'2013-07-01 11:17:51',0),(39,37,1,1629467,'qc1hip5533822','80a14ddb-f9dd-4ad3-b51e-0fd52040bae5',3,2,'2013-07-01 11:17:51',1),(40,40,0,NULL,'qc1hip5529686','qc1hip5529686',NULL,NULL,'2013-07-01 11:17:51',0),(41,40,1,NULL,'qc1hip5529686','qc1hip5529686',NULL,NULL,'2013-07-01 11:17:51',0),(42,40,1,1625186,'qc1hip5529686','59e7eb47-c910-4ca8-9cf0-d829539db05b',4,2,'2013-07-01 11:17:51',1),(43,43,0,NULL,'qc1hip5533821','qc1hip5533821',NULL,NULL,'2013-07-01 11:17:51',0),(44,43,1,NULL,'qc1hip5533821','qc1hip5533821',NULL,NULL,'2013-07-01 11:17:51',0),(45,43,1,1629466,'qc1hip5533821','eacc6a58-e2de-419e-87b7-4f1ce01e8434',3,2,'2013-07-01 11:17:51',1),(46,46,0,NULL,'qc1hip5529685','qc1hip5529685',NULL,NULL,'2013-07-01 11:17:51',0),(47,46,1,NULL,'qc1hip5529685','qc1hip5529685',NULL,NULL,'2013-07-01 11:17:51',0),(48,46,1,1625185,'qc1hip5529685','face6d88-7e90-4215-aa80-fb2c3df5a4ed',5,1,'2013-07-01 11:17:51',1),(49,49,0,NULL,'qc1hip5533833','qc1hip5533833',NULL,NULL,'2013-07-01 11:17:51',0),(50,49,1,NULL,'qc1hip5533833','qc1hip5533833',NULL,NULL,'2013-07-01 11:17:51',0),(51,49,1,1629478,'qc1hip5533833','d95473c2-7e2f-42a1-b6e2-a4beb7c2348b',5,2,'2013-07-01 11:17:51',1),(52,52,0,NULL,'qc1hip5529684','qc1hip5529684',NULL,NULL,'2013-07-01 11:17:51',0),(53,52,1,NULL,'qc1hip5529684','qc1hip5529684',NULL,NULL,'2013-07-01 11:17:51',0),(54,52,1,1625184,'qc1hip5529684','4e9503f1-bd2b-4c59-b09c-3a4060f80901',5,2,'2013-07-01 11:17:51',1),(55,55,0,NULL,'qc1hip5529689','qc1hip5529689',NULL,NULL,'2013-07-01 11:17:51',0),(56,55,1,NULL,'qc1hip5529689','qc1hip5529689',NULL,NULL,'2013-07-01 11:17:51',0),(57,55,1,1625189,'qc1hip5529689','4689a6b5-4dc3-4532-a546-565b41eaa8ba',4,1,'2013-07-01 11:17:51',1),(58,58,0,NULL,'qc1hip5529683','qc1hip5529683',NULL,NULL,'2013-07-01 11:17:52',0),(59,58,1,NULL,'qc1hip5529683','qc1hip5529683',NULL,NULL,'2013-07-01 11:17:52',0),(60,58,1,1625183,'qc1hip5529683','99e5a955-db3c-4010-9b2b-5a68588d0381',5,2,'2013-07-01 11:17:52',1),(61,61,0,NULL,'qc1hip5539494','qc1hip5539494',NULL,NULL,'2013-07-01 11:17:52',0),(62,61,4,NULL,'qc1hip5539494','qc1hip5539494',NULL,NULL,'2013-07-01 11:17:52',0),(63,61,4,1635284,'qc1hip5539494','e3ce4ac6-4d93-4595-8cfa-29bf350d7ac6',6,1,'2013-07-01 11:17:52',1),(64,64,0,NULL,'qc1hip5539493','qc1hip5539493',NULL,NULL,'2013-07-01 11:17:52',0),(65,64,4,NULL,'qc1hip5539493','qc1hip5539493',NULL,NULL,'2013-07-01 11:17:52',0),(66,64,4,1635283,'qc1hip5539493','452d60ca-69d1-4ee4-bf06-ba9f46829537',6,2,'2013-07-01 11:17:52',1),(67,67,0,NULL,'qc1hip5539488','qc1hip5539488',NULL,NULL,'2013-07-01 11:17:52',0),(68,67,4,NULL,'qc1hip5539488','qc1hip5539488',NULL,NULL,'2013-07-01 11:17:52',0),(69,67,4,1635278,'qc1hip5539488','2f2a5ec7-61aa-425f-ba34-07ca8d600129',7,1,'2013-07-01 11:17:52',1),(70,70,0,NULL,'qc1hip5539487','qc1hip5539487',NULL,NULL,'2013-07-01 11:17:52',0),(71,70,4,NULL,'qc1hip5539487','qc1hip5539487',NULL,NULL,'2013-07-01 11:17:52',0),(72,70,4,1635277,'qc1hip5539487','257f8c51-e964-4a4b-aec9-b71a5b10defb',7,2,'2013-07-01 11:17:52',1),(73,73,0,NULL,'qc1hip5543904','qc1hip5543904',NULL,NULL,'2013-07-01 11:17:52',0),(74,73,4,NULL,'qc1hip5543904','qc1hip5543904',NULL,NULL,'2013-07-01 11:17:52',0),(75,73,4,1639801,'qc1hip5543904','7b772883-d419-4623-9428-7b137a84730a',6,2,'2013-07-01 11:17:52',1),(76,76,0,NULL,'qc1hip5539486','qc1hip5539486',NULL,NULL,'2013-07-01 11:17:52',0),(77,76,4,NULL,'qc1hip5539486','qc1hip5539486',NULL,NULL,'2013-07-01 11:17:52',0),(78,76,4,1635276,'qc1hip5539486','e99f1d33-c674-45d4-a96c-f25d883bc9e7',7,2,'2013-07-01 11:17:52',1),(79,79,0,NULL,'qc1hip5543903','qc1hip5543903',NULL,NULL,'2013-07-01 11:17:52',0),(80,79,4,NULL,'qc1hip5543903','qc1hip5543903',NULL,NULL,'2013-07-01 11:17:52',0),(81,79,4,1639800,'qc1hip5543903','7e2c28bf-a861-4268-9f40-5bf442a61158',6,2,'2013-07-01 11:17:52',1),(82,82,0,NULL,'qc1hip5539485','qc1hip5539485',NULL,NULL,'2013-07-01 11:17:52',0),(83,82,4,NULL,'qc1hip5539485','qc1hip5539485',NULL,NULL,'2013-07-01 11:17:52',0),(84,82,4,1635275,'qc1hip5539485','c44d122c-511c-4174-84f2-3333a02df018',7,2,'2013-07-01 11:17:52',1);
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
INSERT INTO `study` VALUES (1,'','2624',NULL,NULL);
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

-- Dump completed on 2013-07-01 10:28:32
