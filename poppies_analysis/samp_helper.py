### Created 6/2/25
## Farhan Hasan (fhasan@stsci.edu)

"""
XPA to SAMP Migration Helper
Convert XPA-based DS9 communication to SAMP using Astropy

This module provides a SAMPHelper class that replaces common XPA operations
with SAMP equivalents for communicating with SAOImage DS9.
"""

import os
import time
import urllib.parse
from pathlib import Path
from astropy.samp import SAMPIntegratedClient
from astropy.io import fits
from astropy.table import Table

import logging

# Set up logging
logger = logging.getLogger(__name__)


class DS9SAMPHelper:
    """Helper class for DS9 SAMP communication."""
    
    def __init__(self, timeout=10):
        """Initialize SAMP client."""
        self.client = SAMPIntegratedClient()
        self.timeout = timeout
        self.connected = False
        self._connect()
    
    def _connect(self):
        """Connect to SAMP hub."""
        try:
            self.client.connect()
            self.connected = True
            logger.info("Connected to SAMP hub successfully")
        except Exception as e:
            logger.warning(f"Failed to connect to SAMP hub: {e}")
            logger.warning("Make sure DS9 is running with SAMP enabled")
            self.connected = False
    
    def _ensure_connected(self):
        """Ensure we're connected, try to reconnect if not."""
        if not self.connected:
            self._connect()
        if not self.connected:
            raise ConnectionError("Cannot connect to SAMP hub. Is DS9 running with SAMP enabled?")
    
    def get_ds9_clients(self):
        """Get list of DS9 clients."""
        self._ensure_connected()
        try:
            all_clients = self.client.get_registered_clients()
            ds9_clients = []
            for client_id in all_clients:
                if client_id == 'hub':
                    continue
                try:
                    metadata = self.client.get_metadata(client_id)
                    # if ('samp.name' in metadata and 
                    #     'ds9' in metadata['samp.name'].lower()):
                    if ('samp.name' in metadata):
                        ds9_clients.append(client_id)
                except:
                    continue
            return ds9_clients
        except Exception as e:
            logger.error(f"Error getting DS9 clients: {e}")
            return []
            
    def send_image_to_ds9(self, fits_path, image_name=None, client_id=None, extension=None):
        """Sends a FITS image to DS9 via SAMP."""
        self._ensure_connected()
        
        fits_path = Path(fits_path).resolve()
        if not fits_path.exists():
            raise FileNotFoundError(f"FITS file not found: {fits_path}")
        
        # Build URL with optional extension specification
        url = fits_path.as_uri()
        if extension is not None:
            url += f"[{extension}]"
        
        params = {
            "url": url,
            "name": image_name or fits_path.name
        }
        
        message = {
            "samp.mtype": "image.load.fits",
            "samp.params": params
        }
        
        try:
            if client_id:
                self.client.notify(client_id, message)
            else:
                self.client.notify_all(message)
            logger.info(f"Sent image {fits_path.name} to DS9")
            
            time.sleep(0.05)

            return True
        
        except Exception as e:
            logger.error(f"Error sending image to DS9: {e}")
            return False
    
    # ### OLDER:
    # def send_image_to_ds9(self, fits_path, image_name=None, client_id=None):
    #     """Sends a FITS image to DS9 via SAMP."""
    #     self._ensure_connected()
        
    #     # print("FH ", self.client.get_registered_clients())
        
    #     # ds9_clients = self.get_ds9_clients()
                
    #     # print("FH testing\n")
    #     # for client_id, meta in ds9_clients.items():
    #     #     print(f"{client_id}: {meta}")
            
    #     fits_path = Path(fits_path).resolve()
    #     if not fits_path.exists():
    #         raise FileNotFoundError(f"FITS file not found: {fits_path}")
        
    #     params = {
    #         "url": fits_path.as_uri(),
    #         "name": image_name or fits_path.name
    #     }
        
    #     message = {
    #         "samp.mtype": "image.load.fits",
    #         "samp.params": params
    #     }
        
    #     try:
    #         if client_id:
    #             self.client.notify(client_id, message)
    #         else:
    #             self.client.notify_all(message)
    #         logger.info(f"Sent image {fits_path.name} to DS9")
    #         return True
    #     except Exception as e:
    #         logger.error(f"Error sending image to DS9: {e}")
    #         return False
    
    def send_command_to_ds9(self, command, client_id=None):
        """Sends a command to DS9 (DS9-specific SAMP extension)."""
        self._ensure_connected()
        
        params = {"cmd": command}
        message = {
            "samp.mtype": "ds9.set",
            "samp.params": params
        }
        
        try:
            if client_id:
                self.client.notify(client_id, message)
            else:
                self.client.notify_all(message)
            logger.debug(f"Sent command to DS9: {command}")

            time.sleep(0.05)

            return True
        
        except Exception as e:
            logger.warning(f"Error sending command '{command}' to DS9: {e}")
            return False
        
    # def identify_ds9_clients_with_labels(base_label='client', delay=0.5):

    def identify_client_labels(self):
        ''' This is the function used for mapping window title to client ID '''

        # client = SAMPIntegratedClient()
        # client.connect()

        self._ensure_connected()

        ds9_clients = {"client_id":[], "name":[]}

        for client_id in self.get_ds9_clients():
            meta = self.client.get_metadata(client_id)
            # print("META: ", client_id, meta)
            if meta.get("samp.name", "").lower().startswith("poppies"):
                name = meta.get("samp.name", "")
                ds9_clients["client_id"].append(client_id)
                ds9_clients["name"].append(name)

    #         print(f"Sent label '{label}' to SAMP client ID: {ds9_id}")
    #         label_map[ds9_id] = label
        # time.sleep(delay)  # Give user time to see it update

        # client.disconnect()
    
        return ds9_clients

    ### ### ### 

    def disconnect(self):
        """Disconnect from SAMP hub."""
        if self.connected:
            try:
                self.client.disconnect()
                self.connected = False
                logger.info("Disconnected from SAMP hub")
            except Exception as e:
                logger.error(f"Error disconnecting: {e}")


def get_ds9_samp_helper():
    """Get or create DS9 SAMP helper instance."""
    global _ds9_samp
    if _ds9_samp is None or not _ds9_samp.connected:
        _ds9_samp = DS9SAMPHelper()
    return _ds9_samp
