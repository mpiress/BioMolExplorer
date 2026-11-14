"""Flet-compatible remembered-device persistence with an in-memory test fallback."""
import flet as ft

KEY='biomolexplorer.auth.remembered_device.v1'


class SessionPersistence:
    def __init__(self,page):
        self.page=getattr(page,'_page',page)
        self.service=None
        self.memory={}
        if isinstance(self.page,ft.Page):
            self.service=ft.SharedPreferences()
            self.page.services.append(self.service)

    async def get(self):
        return await self.service.get(KEY) if self.service else self.memory.get(KEY)

    async def set(self,value):
        if self.service:await self.service.set(KEY,value)
        else:self.memory[KEY]=value

    async def remove(self):
        if self.service:await self.service.remove(KEY)
        else:self.memory.pop(KEY,None)
